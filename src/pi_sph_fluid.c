#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <pthread.h>
#include <time.h>

#include <ssd1306.h>

#define DT 2.5e-4f         // s, time step (stability condition proportional to H/C?)
#define H 0.1300f        // m, particle spacing
#define R (H/1.3f)       // m, initial spacing
#define WIDTH 4.0f       // m, width of domain
#define HEIGHT 2.0f      // m, height of domain
#define RHO_0 1000.0f    // kg/m^3, reference density
#define C 160.0f         // m/s, "numerical" speed of sound (10*max_speed for correct WCSPH)
#define G 9.81f          // m/s^2, gravitational acceleration

#define M (RHO_0*0.83*H*H) // kg, mass of each particle
#define MAX_POSSIBLE_NEIGHBORS 48 // the sum of the first three hexagonal numbers is 22, so this should be enough


// PARTICLES
// SPH particles are implemented as a struct of arrays to suit automatic vectorization

struct particles{
    int count;
    float *x, *y, *rho;
    float *u, *v;
    float *p; // pressure, which is only a function of rho
};

struct particle{ // struct representing a single particle, which should serve as a function input or output
    float x, y, rho;
    float u, v;
    float p; // pressure, which is only a function of rho
};

// Helper function for allocating all the arrays of the particles struct using a malloc-like syntax
struct particles* alloc_particles(int n_particles){
    struct particles *particles = (struct particles*)malloc(sizeof(struct particles));
    
    particles->count = n_particles;
    particles->x = (float*)malloc(n_particles*sizeof(float));
    particles->y = (float*)malloc(n_particles*sizeof(float));
    particles->rho = (float*)malloc(n_particles*sizeof(float));
    particles->u = (float*)malloc(n_particles*sizeof(float));
    particles->v = (float*)malloc(n_particles*sizeof(float));
    particles->p = (float*)malloc(n_particles*sizeof(float));

    return particles;
}

// Helper function for pulling an individual particle out of the struct of arrays
struct particle particle_at(struct particles *particles, int idx){
    return (struct particle){
        .x = particles->x[idx], .y = particles->y[idx], .rho = particles->rho[idx],
        .u = particles->u[idx], .v = particles->v[idx],
        .p = particles->p[idx]
    };
}


// VECTOR RETURN-ELEMENT
// Helper struct for returning two floats from vector-valued functions
typedef struct { float x, y; } float2;


// THE KERNEL AND ITS DERIVATIVE
// Helper functions for calculating the kernel used in SPH to approximate the integral interpolant (in conjunction with 
//  the particle approximation) along with its derivative, using the coordinates of particles a and b as arguments
// Ex: W_ab = W( euclid_dist() / H )    and    grad_a W_ab = ( dW_dq() * dq_dx_a(), dW_dq() * dq_dy_a() )

float euclid_dist(float x_a, float y_a, float x_b, float y_b){
    float x_ab = x_a-x_b, y_ab = y_a-y_b;
    return sqrt(x_ab*x_ab+y_ab*y_ab);
}

float W(float q){
    const float normalizing_factor = 7/(4*M_PI*H*H);
    float tmp_1 = 1-0.5f*q, tmp_2 = 1+2*q;
    return normalizing_factor*pow(tmp_1, 4)*tmp_2; // Wendland C2 kernel, though cubic spline probably also works
}

float dW_dq(float q){
    const float normalizing_factor = 7/(4*M_PI*H*H);
    float tmp = 1-0.5f*q;
    return normalizing_factor*(-5)*q*pow(tmp, 3); // Wendland C2 kernel
}

float dq_dx_a(float x_a, float y_a, float x_b, float y_b){
    return (x_a-x_b)/euclid_dist(x_a, y_a, x_b, y_b)/H;
}

float dq_dy_a(float x_a, float y_a, float x_b, float y_b){
    return (y_a-y_b)/euclid_dist(x_a, y_a, x_b, y_b)/H;
}

float2 grad_a_W_ab(float x_a, float y_a, float x_b, float y_b){
    float q = euclid_dist(x_a, y_a, x_b, y_b)/H;

    float dW_dq_val = dW_dq(q);
    float dq_dx_a_val = dq_dx_a(x_a, y_a, x_b, y_b);
    float dq_dy_a_val = dq_dy_a(x_a, y_a, x_b, y_b);

    return (float2){ .x = dW_dq_val*dq_dx_a_val, .y = dW_dq_val*dq_dy_a_val };
}


// LINKED LIST
// Bare minimum impl of a linked list, along with a couple helper functions, for use in the neighbors search

struct linked_list_element{
    int idx;
    struct linked_list_element *next;
};

struct linked_list{
    struct linked_list_element *head, *tail;
};

void append_element(struct linked_list *list, struct linked_list_element *new_element){
    new_element->next = NULL;

    if(list->head == NULL){
        list->head = new_element;
        list->tail = new_element;
    }
    else{
        list->tail->next = new_element;
        list->tail = new_element;
    }
}


// NEIGHBORS SEARCH
// Old and new implementations of SPH (see Monaghan 2005 and PySPH) recognize that, when the kernel has a "compact 
//  support" (here, W_ab = 0 when euclid_dist() > some number), only a few nearby particles contribute to the particle 
//  approximation of the integral interpolant. This obviously means less operations to do, but implementations 
//  also recognize that finding these nearby particles in the first place could be done more efficiently.
// Here we'll use the linked-list method (mentions in Monaghan 1992 and PySPH). Our kernel's support is a circle of 
//  radius 2*H. Therefore, if a particle falls in some cell of a grid of length 2*H, all possible neighbors are in that 
//  cell and neighboring cells. For each cell, an associated linked list holds the indices of the particles that fall 
//  in it.

struct neighbors_context{
    float x_min, x_max, y_min, y_max; // domain boundaries
    float cell_length;
    int n_cells, m_cells;
    struct linked_list *cells;
    int n_particles;
    struct linked_list_element *cells_elements; // unlike ordinary linked lists, the elements are together alloc'ed once
};

struct neighbors_context *initialize_neighbors_context(int n_particles, float x_min, float x_max, float y_min, 
    float y_max, float cell_length)
{
    struct neighbors_context *context = (struct neighbors_context*)malloc(sizeof(struct neighbors_context));
    
    context->x_min = x_min;
    context->x_max = x_max;
    context->y_min = y_min;
    context->y_max = y_max;
    context->cell_length = cell_length;

    context->n_cells = (int)((y_max-y_min)/cell_length)+1;
    context->m_cells = (int)((x_max-x_min)/cell_length)+1;
    context->cells = (struct linked_list*)malloc(context->n_cells * context->m_cells * sizeof(struct linked_list));
    
    context->n_particles = n_particles;
    context->cells_elements = (struct linked_list_element*)malloc(n_particles*sizeof(struct linked_list_element));
    for(int k = 0; k < n_particles; k++)
        context->cells_elements[k] = (struct linked_list_element){ .idx = k, .next = NULL };
    
    return context;
}

void update_neighbors_context(struct neighbors_context *context, struct particles *particles){
    // reset the linked lists (note that this doesn't orphan the elements)
    for(int ij = 0; ij < context->n_cells * context->m_cells; ij++)
        context->cells[ij] = (struct linked_list){ .head = NULL, .tail = NULL };

    // for each particle, infer the cell it falls in
    for(int k = 0; k < particles->count; k++){
        int i = (int)((particles->y[k] - context->y_min) / context->cell_length);
        int j = (int)((particles->x[k] - context->x_min) / context->cell_length);
        int ij = i * context->m_cells + j;

        append_element(&context->cells[ij], &context->cells_elements[k]);
    }
}

int find_neighbors(int *neighbors_idxs, struct particles *a, struct particles *b, int idx_a, 
    struct neighbors_context *context_b)
{
    int ignore_self_interaction = (a != b); // if a == b (equal ptrs), we need to reject the particle neighboring itself

    // Out of the neighboring cells AND the cell the particle falls in, find the real neighbors
    int neighbors_counter = 0;
    float x_a = a->x[idx_a], y_a = a->y[idx_a];
    int i_a = (int)((y_a - context_b->y_min) / context_b->cell_length),
        j_a = (int)((x_a - context_b->x_min) / context_b->cell_length);
    for(int i = i_a-1; i <= i_a+1; i++){
        for(int j = j_a-1; j <= j_a+1; j++){
            if(i < 0 || i >= context_b->n_cells || j < 0 || j >= context_b->m_cells)
                continue;
            
            int ij = i * context_b->m_cells + j;
            struct linked_list_element *current_element = context_b->cells[ij].head;
            while(current_element != NULL){
                int idx_b = current_element->idx;
                float x_b = b->x[idx_b], y_b = b->y[idx_b];

                if(euclid_dist(x_a, y_a, x_b, y_b) < 2*H && (ignore_self_interaction || idx_a != idx_b)){
                    neighbors_idxs[neighbors_counter] = idx_b;
                    neighbors_counter++;
                }
                
                current_element = current_element->next;
            }
        }
    }

    return neighbors_counter;
}


// SPH APPROXIMATIONS
// Contained implementations of the SPH approximation, integrating all of the above

enum leading_factor { MASS, VOLUME }; // fundamental SPH approx uses volume, but most derived invocations use mass

float sph(float *quantity, struct particles *a, struct particles *b, int idx_a, int *idxs_b, int n_idxs_b,
    enum leading_factor leading_factor)
{
    float x_a = a->x[idx_a], y_a = a->y[idx_a];

    float sph_quantity = 0;
    for(int k = 0; k < n_idxs_b; k++){
        int j = idxs_b[k]; // j is traditionally the index of the neighbor particle

        float x_b = b->x[j], y_b = b->y[j];
        float rho_b = b->rho[j];

        float q = euclid_dist(x_a, y_a, x_b, y_b) / H;
        float W_ab_ij = W(q);
        float leading_factor_j = (leading_factor == MASS)? M : M/rho_b;

        sph_quantity += leading_factor_j * quantity[k] * W_ab_ij;
    }

    return sph_quantity;
}

float2 sph_gradient(float *quantity, struct particles *a, struct particles *b, int idx_a, int *idxs_b, int n_idxs_b,
    enum leading_factor leading_factor)
{
    float x_a = a->x[idx_a], y_a = a->y[idx_a];
    
    float2 grad_quantity = (float2){ .x = 0, .y = 0 };
    for(int k = 0; k < n_idxs_b; k++){
        int j = idxs_b[k];

        float x_b = b->x[j], y_b = b->y[j];
        float rho_b = b->rho[j];

        float2 grad_a_W_ab_ij = grad_a_W_ab(x_a, y_a, x_b, y_b);
        float leading_factor_j = (leading_factor == MASS)? M : M/rho_b;

        grad_quantity.x += leading_factor_j*quantity[k]*grad_a_W_ab_ij.x;
        grad_quantity.y += leading_factor_j*quantity[k]*grad_a_W_ab_ij.y;
    }

    return grad_quantity;
}

float sph_divergence(float *quantity_x, float *quantity_y, struct particles *a, struct particles *b, int idx_a, 
    int *idxs_b, int n_idxs_b, enum leading_factor leading_factor)
{
    float x_a = a->x[idx_a], y_a = a->y[idx_a];

    float div_quantity = 0;
    for(int k = 0; k < n_idxs_b; k++){
        int j = idxs_b[k];

        float x_b = b->x[j], y_b = b->y[j];
        float rho_b = b->rho[j];

        float2 grad_a_W_ab_ij = grad_a_W_ab(x_a, y_a, x_b, y_b);
        float leading_factor_j = (leading_factor == MASS)? M : M/rho_b;

        float quantity_dot_grad = quantity_x[k]*grad_a_W_ab_ij.x + quantity_y[k]*grad_a_W_ab_ij.y;

        div_quantity += leading_factor_j*quantity_dot_grad;
    }

    return div_quantity;
}


// MAIN FUNCTIONS
// These functions are responsible for principal parts of the fluid simulation. Exact implementation of Monaghan 1992 
//  to the best of my ability with the exception of leaving artificial viscosity on all the time

// To start, arrange the SPH particles as a circle in the middle of the domain
int in_initial_shape(float x, float y){
    return euclid_dist(x, y, WIDTH/2, HEIGHT/2) < 0.70;
}

void calculate_density(struct particles *fluid, struct particles *boundary, struct neighbors_context *ctx_fluid, 
    struct neighbors_context *ctx_boundary)
{
    static float ones[MAX_POSSIBLE_NEIGHBORS] = {1, 1, 1, 1, 1, 1, 1, 1, 
        1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 
        1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
    int neighbors_idxs[MAX_POSSIBLE_NEIGHBORS], n_neighbors;

    #pragma omp for
    for(int i = 0; i < fluid->count; i++){
        float rho_i = 1.293; // init with a small arbitrary val (chose air density) to avoid div-by-zero
        
        // fluid contribution to density of fluid
        n_neighbors = find_neighbors(neighbors_idxs, fluid, fluid, i, ctx_fluid);
        rho_i += sph(ones, fluid, fluid, i, neighbors_idxs, n_neighbors, MASS);

        // boundary contribution to density of fluid
        n_neighbors = find_neighbors(neighbors_idxs, fluid, boundary, i, ctx_boundary);
        rho_i += sph(ones, fluid, boundary, i, neighbors_idxs, n_neighbors, MASS);

        fluid->rho[i] = rho_i;
    }

    #pragma omp for
    for(int b = 0; b < boundary->count; b++){
        float rho_b = 1.293;

        // fluid contribution to density of boundary
        n_neighbors = find_neighbors(neighbors_idxs, boundary, fluid, b, ctx_fluid);
        rho_b += sph(ones, boundary, fluid, b, neighbors_idxs, n_neighbors, MASS);

        boundary->rho[b] = rho_b;
    }
}

// "slighly compressible" SPH, or weakly-compressible SPH (WCSPH) expresses pressure as an explicit function of density,
//  not an implicit one needing an iterative solver. Sometimes doesn't use the true speed of sound in a fluid, and
//  rather uses some number high enough that the density doesn't vary by too many percents. That's also how I did it
//  here because I saw that using the true speed caused instability. See Monaghan 1992 or Monaghan 2005.
void calculate_particle_pressure(struct particles *particles){
    #pragma omp for
    for(int i = 0; i < particles->count; i++){
        const float B = C*C*RHO_0/7;
        float rho_ratio = particles->rho[i]/RHO_0;
        float pressure_i = B*(pow(rho_ratio, 7)-1);
        particles->p[i] = (pressure_i > 0)? pressure_i : 0;
    }
}

void add_pressure_acceleration(float *du_dt_fluid, float *dv_dt_fluid, struct particles *fluid, 
    struct particles *boundary, struct neighbors_context *ctx_fluid, struct neighbors_context *ctx_boundary)
{
    int neighbors_idxs[MAX_POSSIBLE_NEIGHBORS], n_neighbors;
    float pressure_i[MAX_POSSIBLE_NEIGHBORS];
    float2 pressure_grad_i;
    
    #pragma omp for
    for(int idx_a = 0; idx_a < fluid->count; idx_a++){
        struct particle a_i = particle_at(fluid, idx_a);       
         
        float du_dt_fluid_i = 0, dv_dt_fluid_i = 0;

        // fluid contribution to pressure acceleration of fluid
        n_neighbors = find_neighbors(neighbors_idxs, fluid, fluid, idx_a, ctx_fluid);

        // compute parts of the momentum-conserving pressure from neighbors
        for(int k = 0; k < n_neighbors; k++){
            int idx_b = neighbors_idxs[k];
            struct particle b_j = particle_at(fluid, idx_b);

            pressure_i[k] = -( a_i.p/(a_i.rho*a_i.rho) + b_j.p/(b_j.rho*b_j.rho) );
        }

        // compute the acceleration due to pressure using the SPH gradient
        pressure_grad_i = sph_gradient(pressure_i, fluid, fluid, idx_a, neighbors_idxs, n_neighbors, MASS);
        du_dt_fluid_i += pressure_grad_i.x;
        dv_dt_fluid_i += pressure_grad_i.y;

        // boundary contribution to pressure acceleration of fluid
        n_neighbors = find_neighbors(neighbors_idxs, fluid, boundary, idx_a, ctx_boundary);

        for(int k = 0; k < n_neighbors; k++){
            int idx_b = neighbors_idxs[k];
            struct particle b_j = particle_at(boundary, idx_b);

            pressure_i[k] = -( a_i.p/(a_i.rho*a_i.rho) + b_j.p/(b_j.rho*b_j.rho) );
        }

        pressure_grad_i = sph_gradient(pressure_i, fluid, boundary, idx_a, neighbors_idxs, n_neighbors, MASS);
        du_dt_fluid_i += pressure_grad_i.x;
        dv_dt_fluid_i += pressure_grad_i.y;

        du_dt_fluid[idx_a] += du_dt_fluid_i;
        dv_dt_fluid[idx_a] += dv_dt_fluid_i;
    }
}

void add_viscosity_acceleration(float *du_dt_fluid, float *dv_dt_fluid, struct particles *fluid, 
    struct particles *boundary, struct neighbors_context *ctx_fluid, struct neighbors_context *ctx_boundary)
{
    int neighbors_idxs[MAX_POSSIBLE_NEIGHBORS], n_neighbors;
    float viscosity_i[MAX_POSSIBLE_NEIGHBORS];
    float2 viscosity_grad_i;
    
    #pragma omp for
    for(int idx_a = 0; idx_a < fluid->count; idx_a++){
        struct particle a_i = particle_at(fluid, idx_a);       
         
        float du_dt_fluid_i = 0, dv_dt_fluid_i = 0;

        // fluid contribution to pressure acceleration of fluid
        n_neighbors = find_neighbors(neighbors_idxs, fluid, fluid, idx_a, ctx_fluid);

        // compute parts of the momentum-conserving pressure from neighbors
        for(int k = 0; k < n_neighbors; k++){
            int idx_b = neighbors_idxs[k];
            struct particle b_j = particle_at(fluid, idx_b);

            float u_ab = a_i.u-b_j.u, v_ab = a_i.v-b_j.v;
            float x_ab = a_i.x-b_j.x, y_ab = a_i.y-b_j.y;
            float mean_rho = (a_i.rho+b_j.rho)/2;

            float xy_dot_uv = x_ab*u_ab+y_ab*v_ab;
            float xy_dot_xy = x_ab*x_ab+y_ab*y_ab;
            float mu_ab = H*xy_dot_uv/(xy_dot_xy+0.01*H*H);

            viscosity_i[k] = 0.01*C*mu_ab/mean_rho;
        }

        // compute the acceleration due to pressure using the SPH gradient
        viscosity_grad_i = sph_gradient(viscosity_i, fluid, fluid, idx_a, neighbors_idxs, n_neighbors, MASS);
        du_dt_fluid_i += viscosity_grad_i.x;
        dv_dt_fluid_i += viscosity_grad_i.y;

        // boundary contribution to pressure acceleration of fluid
        n_neighbors = find_neighbors(neighbors_idxs, fluid, boundary, idx_a, ctx_boundary);

        for(int k = 0; k < n_neighbors; k++){
            int idx_b = neighbors_idxs[k];
            struct particle b_j = particle_at(boundary, idx_b);

            float u_ab = a_i.u-b_j.u, v_ab = a_i.v-b_j.v;
            float x_ab = a_i.x-b_j.x, y_ab = a_i.y-b_j.y;
            float mean_rho = (a_i.rho+b_j.rho)/2;

            float xy_dot_uv = x_ab*u_ab+y_ab*v_ab;
            float xy_dot_xy = x_ab*x_ab+y_ab*y_ab;
            float mu_ab = H*xy_dot_uv/(xy_dot_xy+0.01*H*H);

            viscosity_i[k] = 0.01*C*mu_ab/mean_rho;
        }

        viscosity_grad_i = sph_gradient(viscosity_i, fluid, boundary, idx_a, neighbors_idxs, n_neighbors, MASS);
        du_dt_fluid_i += viscosity_grad_i.x;
        dv_dt_fluid_i += viscosity_grad_i.y;

        du_dt_fluid[idx_a] += du_dt_fluid_i;
        dv_dt_fluid[idx_a] += dv_dt_fluid_i;
    }
}

void add_compression_effect(float *drho_dt, struct particles *a, struct particles *b, 
    struct neighbors_context *context_b)
{
    #pragma omp for
    for(int idx_a = 0; idx_a < a->count; idx_a++){
        struct particle a_i = particle_at(a, idx_a);

        int neighbors_idxs[MAX_POSSIBLE_NEIGHBORS];
        int n_neighbors = find_neighbors(neighbors_idxs, a, b, idx_a, context_b);

        // compute parts of the popular formulation of the continuity equation from neighbors
        float u_ab[MAX_POSSIBLE_NEIGHBORS], v_ab[MAX_POSSIBLE_NEIGHBORS];
        for(int k = 0; k < n_neighbors; k++){
            int idx_b = neighbors_idxs[k];
            struct particle b_j = particle_at(b, idx_b);

            u_ab[k] = a_i.u-b_j.u;
            v_ab[k] = a_i.v-b_j.v;
        }

        // compute the change in density using the SPH divergence
        float drho_dt_i = sph_divergence(u_ab, v_ab, a, b, idx_a, neighbors_idxs, n_neighbors, MASS);
        drho_dt[idx_a] += drho_dt_i;
    }
}


// METABALLS
// See Wikipedia article and original Blinn 1982 paper on metaballs. This function follows from the original derivation 
//  of a surface from points, but it's currently not properly parameterized or optimized. On the other hand, it manages 
//  to hook off of the existing neighbor-finding code and kernel function.

void draw_metaballs(unsigned char *draw_buffer, struct particles *pixel_pseudoparticles, struct particles *fluid, 
    struct neighbors_context *ctx_fluid)
{
    #pragma omp for collapse(2)
    for(int i = 0; i < 64; i++){
        for(int j = 0; j < 128; j++){
            float pixel_x = pixel_pseudoparticles->x[i*128+j],
                    pixel_y = pixel_pseudoparticles->y[i*128+j];
            int ij = i*128+j;

            int contributing_particles[MAX_POSSIBLE_NEIGHBORS];
            int n_contributing_particles = find_neighbors(contributing_particles, pixel_pseudoparticles, fluid, ij, ctx_fluid);

            float metaball_condition = 0;
            for(int k = 0; k < n_contributing_particles; k++){
                int particle_index = contributing_particles[k];
                float particle_x = fluid->x[particle_index],
                        particle_y = fluid->y[particle_index];
                
                float distance = euclid_dist(particle_x, particle_y, pixel_x, pixel_y);
                float q = distance/H;

                float unnormalizing_factor = (4*M_PI*H*H)/7;
                float new_normalizing_factor = 0.75;
                metaball_condition += new_normalizing_factor*unnormalizing_factor*W(q);

                if(metaball_condition >= 1) break;
            }

            #pragma omp critical
            if(metaball_condition >= 1) draw_buffer[i/8*128+j] |= (1 << (i%8));
            else draw_buffer[i/8*128+j] &= ~(1 << (i%8));
        }
    }
}


// AUXILIARY ROUTINES
// Interacts with the simulation via MPU6050 accelerometer and SSD1306 OLED display. Runs on separate pthreads.

int read_file_as_integer(const char *filepath){
    FILE *file = fopen(filepath, "r");
    if(file == NULL){
        printf("Error opening file %s\n", filepath);
        exit(1);
    }

    int value;
    fscanf(file, "%d", &value);
    fclose(file);
    return value;
}

// For now, reads and interprets values generated by the driver integreted in the Linux kernel (TODO: userspace driver?)
// (init driver with "echo mpu6050 0x68 > /sys/bus/i2c/devices/i2c-1/new_device")
void get_gravity(float *g_x, float *g_y){
    #ifdef MPU6050
    int accel_x = read_file_as_integer("/sys/bus/iio/devices/iio:device0/in_accel_x_raw"),
        accel_y = read_file_as_integer("/sys/bus/iio/devices/iio:device0/in_accel_y_raw"),
        accel_z = read_file_as_integer("/sys/bus/iio/devices/iio:device0/in_accel_z_raw");
    float accel_magnitude = sqrt(accel_x*accel_x+accel_y*accel_y+accel_z*accel_z);
    float normed_accel_x = accel_x/accel_magnitude,
          normed_accel_y = accel_y/accel_magnitude;
    *g_x = normed_accel_y*G;
    *g_y = -normed_accel_x*G;
    #else
    *g_x = 0;
    *g_y = -G;
    #endif
}

void* get_gravity_routine(void *arg){
    float *g = (float*)arg;
    
    struct timespec last, now;
    clock_gettime(CLOCK_MONOTONIC, &last);
    clock_gettime(CLOCK_MONOTONIC, &now);
    while(1){
        while(now.tv_nsec-last.tv_nsec < 1000000000/20){
            usleep(1000);
            clock_gettime(CLOCK_MONOTONIC, &now);
        }

        get_gravity(&g[0], &g[1]);
    }
}

void* display_routine(void *arg){
    unsigned char *display_buffer = (unsigned char*)arg;

    struct timespec last, now;
    clock_gettime(CLOCK_MONOTONIC, &last);
    clock_gettime(CLOCK_MONOTONIC, &now);
    while(1){
        while(now.tv_nsec-last.tv_nsec < 1000000000/60){
            usleep(1000);
            clock_gettime(CLOCK_MONOTONIC, &now);
        }

        ssd1306_drawBufferFast(0, 0, 128, 64, display_buffer);
    }
}


// MAIN

int main(){
    int particle_counter; // in case the number of particles is determined numerically, this is used to count them


    // initialize display
    ssd1306_128x64_i2c_init();
    ssd1306_clearScreen();
    
    
    // count the number of fluid particles we need
    int n_fluid = 0;
    for(float x_0 = 0; x_0 < WIDTH; x_0 += R)
        for(float y_0 = 0; y_0 < HEIGHT; y_0 += R)
            if(in_initial_shape(x_0, y_0)) n_fluid++;
    
    // alloc fluid and derivatives
    struct particles *fluid, *fluid_pred;
    fluid = alloc_particles(n_fluid);
    fluid_pred = alloc_particles(n_fluid);
    float *du_dt_pred = (float*)malloc(n_fluid*sizeof(float)), // momentum and continuity results for predictor step
          *dv_dt_pred = (float*)malloc(n_fluid*sizeof(float));
    float *du_dt_corr = (float*)malloc(n_fluid*sizeof(float)), // momentum and continuity results for corrector step
          *dv_dt_corr = (float*)malloc(n_fluid*sizeof(float));
    
    // initialize fluid particles
    particle_counter = 0;
    for(float x_0 = 0; x_0 < WIDTH; x_0 += R){
        for(float y_0 = 0; y_0 < HEIGHT; y_0 += R){
            if(in_initial_shape(x_0, y_0)){
                fluid->x[particle_counter] = x_0;
                fluid->y[particle_counter] = y_0;
                fluid->u[particle_counter] = 0;
                fluid->v[particle_counter] = 0;
                fluid->rho[particle_counter] = RHO_0;

                particle_counter++;
            }
        }
    }


    // count the number of boundary particles we need
    int n_boundary = 0;
    for(float x_0 = 0; x_0 < WIDTH; x_0 += R/3) n_boundary += 4; // two particles per side (top and bottom)
    for(float y_0 = 0; y_0 < HEIGHT; y_0 += R/3) n_boundary += 4;

    // alloc boundary particles and derivatives
    struct particles *boundary, *boundary_pred; 
    boundary = alloc_particles(n_boundary);
    boundary_pred = alloc_particles(n_boundary);

    // initialize boundary particles velocity and density (velocity will never change)
    for(int i = 0; i < n_boundary; i++){
        boundary->u[i] = 0;
        boundary->v[i] = 0;
        boundary->rho[i] = RHO_0;
    }

    // initialize boundary particles locations (will never change)
    particle_counter = 0;
    for(float x_0 = 0; x_0 < WIDTH; x_0 += R/3){
        boundary->x[particle_counter] = x_0;
        boundary->y[particle_counter] = 0;
        boundary->x[particle_counter+1] = x_0;
        boundary->y[particle_counter+1] = 0-R;
        boundary->x[particle_counter+2] = x_0;
        boundary->y[particle_counter+2] = HEIGHT;
        boundary->x[particle_counter+3] = x_0;
        boundary->y[particle_counter+3] = HEIGHT+R;
        particle_counter += 4;
    }
    for(float y_0 = 0; y_0 < HEIGHT; y_0 += R/3){
        boundary->x[particle_counter] = 0;
        boundary->y[particle_counter] = y_0;
        boundary->x[particle_counter+1] = 0-R;
        boundary->y[particle_counter+1] = y_0;
        boundary->x[particle_counter+2] = WIDTH;
        boundary->y[particle_counter+2] = y_0;
        boundary->x[particle_counter+3] = WIDTH+R;
        boundary->y[particle_counter+3] = y_0;
        particle_counter += 4;
    }

    // boundary_pred x, y, u, and v are the same as boundary but will never be 
    // updated, so we need to initialize them now
    for(int i = 0; i < n_boundary; i++){
        boundary_pred->x[i] = boundary->x[i];
        boundary_pred->y[i] = boundary->y[i];
        boundary_pred->u[i] = boundary->u[i];
        boundary_pred->v[i] = boundary->v[i];
    }


    printf("n_fluid = %d\n", n_fluid);
    printf("n_boundary = %d\n", n_boundary);


    // initialize neighbors search context
    float x_min = 0-R, x_max = WIDTH+R, y_min = 0-R, y_max = HEIGHT+R;
    struct neighbors_context *ctx_fluid, *ctx_boundary;
    ctx_fluid = initialize_neighbors_context(fluid->count, x_min, x_max, y_min, y_max, 2*H);
    ctx_boundary = initialize_neighbors_context(boundary->count, x_min, x_max, y_min, y_max, 2*H);
    update_neighbors_context(ctx_boundary, boundary); // this never needs to be called again, so we'll call it now


    struct timespec now; // initialize the time-keeping with the current time
    clock_gettime(CLOCK_MONOTONIC, &now);


    // set up before metaballs drawing and launch the display thread
    pthread_t display_thread;
    unsigned char *draw_buffer = (unsigned char*)calloc(1024, 1);
    pthread_create(&display_thread, NULL, display_routine, draw_buffer);

    // in leiu of defining a new function for finding the contributing particles to the metaballs condition, we'll just 
    //   reuse the neighbors search function (called with ctx_fluid as the argument)
    // to do so, we define pseudoparticles at the pixel centers (not unlike how we do treat the boundary)
    struct particles *pixel_pseudoparticles = alloc_particles(64*128);
    for(int i = 0; i < 64; i++){
        for(int j = 0; j < 128; j++){
            float pixel_x = (j+0.5)*WIDTH/128, pixel_y = (64-(i+0.5))*HEIGHT/64;
            pixel_pseudoparticles->x[i*128+j] = pixel_x;
            pixel_pseudoparticles->y[i*128+j] = pixel_y;
        }
    }

    struct timespec last_drew = now; // initialize the last-drew time to now


    // initialize gravity and the gravity-reading thread
    float g[2];
    get_gravity(&g[0], &g[1]);
    pthread_t gravity_thread;
    pthread_create(&gravity_thread, NULL, get_gravity_routine, g);


    // initialize statistics reporting and buffer drawing
    float worst_avg_rho_error_pct = 0;
    float t = 0, last_t = 0;
    struct timespec last_reported = now;
    clock_gettime(CLOCK_MONOTONIC, &last_reported);


    // main loop
    #pragma omp parallel num_threads(4)
    while(1){
        #pragma omp single
        {
            // predictor step: initialize the derivative accumulators
            for(int i = 0; i < n_fluid; i++){
                du_dt_pred[i] = g[0];
                dv_dt_pred[i] = g[1];
            }

            // predictor step: update the neighbors search context
            update_neighbors_context(ctx_fluid, fluid);
        }

        // predictor step: calculate pressure and take the sum of contributions to the derivatives from the neighbors
        calculate_density(fluid, boundary, ctx_fluid, ctx_boundary);
        calculate_particle_pressure(fluid);
        calculate_particle_pressure(boundary);
        add_pressure_acceleration(du_dt_pred, dv_dt_pred, fluid, boundary, ctx_fluid, ctx_boundary);
        add_viscosity_acceleration(du_dt_pred, dv_dt_pred, fluid, boundary, ctx_fluid, ctx_boundary);

        #pragma omp single
        {
            // predictor step: get what the particles would be like if we used forward Euler
            for(int i = 0; i < n_fluid; i++){
                fluid_pred->x[i] = fluid->x[i] + fluid->u[i]*DT;
                fluid_pred->y[i] = fluid->y[i] + fluid->v[i]*DT;
                fluid_pred->u[i] = fluid->u[i] + du_dt_pred[i]*DT;
                fluid_pred->v[i] = fluid->v[i] + dv_dt_pred[i]*DT;
            }

            // corrector step: initialize the derivative accumulators
            for(int i = 0; i < n_fluid; i++){
                du_dt_corr[i] = g[0];
                dv_dt_corr[i] = g[1];
            }

            // corrector step: update the neighbors context using the predictor positions
            update_neighbors_context(ctx_fluid, fluid_pred);
        }

        // corrector step: calculate pressure and take the sum of contributions to the derivatives from the neighbors
        calculate_density(fluid_pred, boundary_pred, ctx_fluid, ctx_boundary);
        calculate_particle_pressure(fluid_pred);
        calculate_particle_pressure(boundary_pred);
        add_pressure_acceleration(du_dt_corr, dv_dt_corr, fluid_pred, boundary_pred, ctx_fluid, ctx_boundary);
        add_viscosity_acceleration(du_dt_corr, dv_dt_corr, fluid_pred, boundary_pred, ctx_fluid, ctx_boundary);

        #pragma omp single
        {
            // corrector step: step forward using the midpoint between the predictor and corrector derivatives
            for(int i = 0; i < n_fluid; i++){
                fluid->x[i] += 0.5f*(fluid_pred->u[i] + fluid->u[i])*DT;
                fluid->y[i] += 0.5f*(fluid_pred->v[i] + fluid->v[i])*DT;
                fluid->u[i] += 0.5f*(du_dt_pred[i] + du_dt_corr[i])*DT;
                fluid->v[i] += 0.5f*(dv_dt_pred[i] + dv_dt_corr[i])*DT;
            }


            clock_gettime(CLOCK_MONOTONIC, &now); // take a single timestamp for the below real-time operations
        }


        // draw fluid using metaballs
        float elapsed = (now.tv_sec-last_drew.tv_sec) + (now.tv_nsec-last_drew.tv_nsec)/1e9;
        if(elapsed > 1.0/60){
            draw_metaballs(draw_buffer, pixel_pseudoparticles, fluid, ctx_fluid);

            #pragma omp single
            last_drew = now;
        }


        #pragma omp single
        {    
            // record the average rho ratio in a single frame as avg_rho_ratio
            // compare avg_rho_ratio to worst_avg_rho_error_pct, which is the worst average rho ratio out of ALL frames
            // this is a critical statistic that shows to what degree the incompressibility constraint is being violated
            float temp = 0, avg_rho_error, avg_rho_error_pct;
            for(int i = 0; i < fluid->count; i++) temp += fluid->rho[i];
            avg_rho_error = temp / fluid->count - RHO_0;
            avg_rho_error_pct = avg_rho_error/RHO_0*100;
            if(avg_rho_error_pct > worst_avg_rho_error_pct) worst_avg_rho_error_pct = avg_rho_error;

            
            // report frame rate and other statistics
            t += DT;
            if(t-last_t > 0.1){
                float elapsed = (now.tv_sec-last_reported.tv_sec) + (now.tv_nsec-last_reported.tv_nsec)/1e9;
                int tps = ((t-last_t)/DT)/elapsed;

                printf("t=%.2f, ", t);
                printf("ticks/s = %d, ", tps);
                printf("worst avg rho error = %.3f\%%, ", worst_avg_rho_error_pct);
                printf("current avg rho error = %.3f\%%, ", avg_rho_error_pct);
                printf("\n");

                last_t = t;
                last_reported = now;
            }
        }
    }
}