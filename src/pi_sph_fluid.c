#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <pthread.h>
#include <time.h>

#include <ssd1306.h>

#define R 0.1000f        // m, initial spacing (real ticks/s is O(R^3), but DT is O(R), so realtime implies intersect)
#define H (R*1.3f)       // m, smoothing length
#define WIDTH 4.0f       // m, width of domain
#define HEIGHT 2.0f      // m, height of domain
#define RHO_0 1000.0f    // kg/m^3, reference density
#define C 160.0f         // m/s, "numerical" speed of sound (10*max_speed for correct WCSPH)
#define G 9.81f          // m/s^2, gravitational acceleration

#define DT (0.4f*H/C)    // s, time step (CFL condition? DT <= 0.4*H/C)
#define V (0.83*H*H)     // m^3, volume of each fluid particle
#define MAX_POSSIBLE_NEIGHBORS 48 // the sum of the first three hexagonal numbers is 22, so this should be enough


// PARTICLES

struct particle{ // struct representing a single particle
    float x, y, u, v; // intrinsic values
    float m; // mass is RHO_0*V for fluid particles, but some calculated constant for boundary particles
    float rho; // rho is derived using SPH
    float p; // pressure is derived from some incompressible-enforcing routine (here being WCSPH)
};


// VECTOR RETURN-ELEMENT
// Helper struct for returning two floats from vector-valued functions
typedef struct { float x, y; } float2;


// THE KERNEL AND ITS DERIVATIVE
// Helper functions for calculating the kernel used in SPH to approximate the integral interpolant (in conjunction with 
//  the particle approximation) along with its derivative, using the coordinates of particles a and b as arguments
// Ex: W_ab = W( euclid_dist() / H )    and    grad_a W_ab = ( dW_dq() * dq_dx_a(), dW_dq() * dq_dy_a() )

float euclid_dist(float x_i, float y_i, float x_j, float y_j){
    float x_ij = x_i-x_j, y_ij = y_i-y_j;
    return sqrt(x_ij*x_ij+y_ij*y_ij);
}

float W(float q){
    const float normalizing_factor = 7/(4*M_PI*H*H);
    float tmp_1 = 1-0.5f*q, tmp_2 = 1+2*q;
    return normalizing_factor*pow(tmp_1, 4)*tmp_2; // Wendland C2 kernel, though cubic spline probably also works
}

float2 grad_a_W_ab(float x_i, float y_i, float x_j, float y_j){
    const float normalizing_factor = 7/(4*M_PI*H*H);
    float q = euclid_dist(x_i, y_i, x_j, y_j)/H;
    float tmp = 1-0.5f*q;
    float dW_dq_ij = normalizing_factor*(-5)*q*pow(tmp, 3); // Wendland C2 kernel
    
    float dq_dx_a_ij = (x_i-x_j)/euclid_dist(x_i, y_i, x_j, y_j)/H;
    float dq_dy_a_ij = (y_i-y_j)/euclid_dist(x_i, y_i, x_j, y_j)/H;

    return (float2){ .x = dW_dq_ij*dq_dx_a_ij, .y = dW_dq_ij*dq_dy_a_ij };
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
// Here we'll use the linked-list method (mentions in Monaghan 1994 and PySPH). Our kernel's support is a circle of 
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

struct neighbors_context *alloc_neighbors_context(int n_particles, float x_min, float x_max, float y_min, 
    float y_max, float cell_length)
{
    struct neighbors_context *ctx = (struct neighbors_context*)malloc(sizeof(struct neighbors_context));
    
    ctx->x_min = x_min;
    ctx->x_max = x_max;
    ctx->y_min = y_min;
    ctx->y_max = y_max;
    ctx->cell_length = cell_length;

    ctx->n_cells = (int)((y_max-y_min)/cell_length)+1;
    ctx->m_cells = (int)((x_max-x_min)/cell_length)+1;
    ctx->cells = (struct linked_list*)malloc(ctx->n_cells * ctx->m_cells * sizeof(struct linked_list));
    
    ctx->n_particles = n_particles;
    ctx->cells_elements = (struct linked_list_element*)malloc(n_particles*sizeof(struct linked_list_element));
    for(int i = 0; i < n_particles; i++)
        ctx->cells_elements[i] = (struct linked_list_element){ .idx = i, .next = NULL };
    
    return ctx;
}

void update_neighbors_context(struct neighbors_context *ctx, struct particle *particles){
    // reset the linked lists (note that this doesn't orphan the elements)
    for(int ij_cell = 0; ij_cell < ctx->n_cells * ctx->m_cells; ij_cell++)
        ctx->cells[ij_cell] = (struct linked_list){ .head = NULL, .tail = NULL };

    // for each particle, infer the cell it falls in
    for(int i = 0; i < ctx->n_particles; i++){
        int i_cell = (int)((particles[i].y - ctx->y_min) / ctx->cell_length);
        int j_cell = (int)((particles[i].x - ctx->x_min) / ctx->cell_length);
        int ij_cell = i_cell * ctx->m_cells + j_cell;

        append_element(&ctx->cells[ij_cell], &ctx->cells_elements[i]);
    }
}

int find_neighbors(int *j_neighbors, struct particle *particles_a, struct particle *particles_b, int i, 
    struct neighbors_context *ctx_b)
{
    // if particles_a == particles_b (equal ptrs), we need to reject the particle neighboring itself
    int ignore_self_interaction = (particles_a != particles_b);

    // Out of the neighboring cells AND the cell the particle falls in, find the real neighbors
    int neighbors_counter = 0;
    int i_cell_center = (int)((particles_a[i].y - ctx_b->y_min) / ctx_b->cell_length),
        j_cell_center = (int)((particles_a[i].x - ctx_b->x_min) / ctx_b->cell_length);
    for(int i_cell = i_cell_center-1; i_cell <= i_cell_center+1; i_cell++){
        for(int j_cell = j_cell_center-1; j_cell <= j_cell_center+1; j_cell++){
            if(i_cell < 0 || i_cell >= ctx_b->n_cells || j_cell < 0 || j_cell >= ctx_b->m_cells)
                continue;
            int ij_cell = i_cell * ctx_b->m_cells + j_cell;

            struct linked_list_element *current_element = ctx_b->cells[ij_cell].head;
            while(current_element != NULL){
                int j = current_element->idx;
                float distance = euclid_dist(particles_a[i].x, particles_a[i].y, particles_b[j].x, particles_b[j].y);

                if(distance < 2*H && (ignore_self_interaction || i != j)){
                    j_neighbors[neighbors_counter] = j;
                    neighbors_counter++;
                }
                
                current_element = current_element->next;
            }
        }
    }

    return neighbors_counter;
}

// neighbors should be organized as a struct of arrays to suit automatic vectorization, but this can otherwise be 
//  thought of as an abstract object containing a copy of the particle data from each neighbor
struct particle_neighbors{
    int count;
    float x[MAX_POSSIBLE_NEIGHBORS], y[MAX_POSSIBLE_NEIGHBORS], u[MAX_POSSIBLE_NEIGHBORS], v[MAX_POSSIBLE_NEIGHBORS];
    float m[MAX_POSSIBLE_NEIGHBORS];
    float rho[MAX_POSSIBLE_NEIGHBORS];
    float p[MAX_POSSIBLE_NEIGHBORS];
};

// to do so, we'll combine the fetch step with a transpose from AoS to SoA
void read_neighbors(struct particle *particles, int *j_neighbors, int n_neighbors, 
    struct particle_neighbors *neighbors)
{
    neighbors->count = n_neighbors;

    for(int k = 0; k < n_neighbors; k++){
        int j = j_neighbors[k];

        neighbors->x[k] = particles[j].x;
        neighbors->y[k] = particles[j].y;
        neighbors->u[k] = particles[j].u;
        neighbors->v[k] = particles[j].v;
        neighbors->m[k] = particles[j].m;
        neighbors->rho[k] = particles[j].rho;
        neighbors->p[k] = particles[j].p;
    }
}

// Helper function for pulling an individual particle out of the struct of arrays
struct particle particle_at(struct particle_neighbors *particles, int i){
    return (struct particle){
        .x = particles->x[i], .y = particles->y[i], .u = particles->u[i], .v = particles->v[i],
        .m = particles->m[i],
        .rho = particles->rho[i],
        .p = particles->p[i]
    };
}


// SPH APPROXIMATIONS
// Contained implementations of the SPH approximation, integrating all of the above

enum leading_factor { MASS, VOLUME }; // fundamental SPH approx uses volume, but most derived invocations use mass

float sph(float *quantity, struct particle particle_i, struct particle_neighbors *particle_i_neighbors, 
    enum leading_factor leading_factor)
{
    float sph_quantity = 0;
    for(int k = 0; k < particle_i_neighbors->count; k++){
        struct particle neighbor_j = particle_at(particle_i_neighbors, k);

        float W_ij = W(euclid_dist(particle_i.x, particle_i.y, neighbor_j.x, neighbor_j.y) / H);
        float leading_factor_j = (leading_factor == MASS)? neighbor_j.m : neighbor_j.m/neighbor_j.rho;

        sph_quantity += leading_factor_j * quantity[k] * W_ij;
    }

    return sph_quantity;
}

float2 sph_gradient(float *quantity, struct particle particle_i, struct particle_neighbors *particle_i_neighbors, 
    enum leading_factor leading_factor)
{
    float2 grad_quantity = (float2){ .x = 0, .y = 0 };
    for(int k = 0; k < particle_i_neighbors->count; k++){
        struct particle neighbor_j = particle_at(particle_i_neighbors, k);

        float2 grad_i_W_ij = grad_a_W_ab(particle_i.x, particle_i.y, neighbor_j.x, neighbor_j.y);
        float leading_factor_j = (leading_factor == MASS)? neighbor_j.m : neighbor_j.m/neighbor_j.rho;

        grad_quantity.x += leading_factor_j*quantity[k]*grad_i_W_ij.x;
        grad_quantity.y += leading_factor_j*quantity[k]*grad_i_W_ij.y;
    }

    return grad_quantity;
}


// MAIN FUNCTIONS
// These functions are responsible for principal parts of the fluid simulation. Exact implementation of Monaghan 1994 
//  to the best of my ability with the exception of leaving artificial viscosity on all the time

// To start, arrange the SPH particles as a circle in the middle of the domain
int in_initial_shape(float x, float y){
    return euclid_dist(x, y, WIDTH/2, HEIGHT/2) < 0.70;
}

void calculate_boundary_pseudomass(struct particle *boundary, struct neighbors_context *ctx_boundary){
    int j_neighbors[MAX_POSSIBLE_NEIGHBORS], n_neighbors;
    struct particle_neighbors neighbors;

    #pragma omp for
    for(int i = 0; i < ctx_boundary->n_particles; i++){
        n_neighbors = find_neighbors(j_neighbors, boundary, boundary, i, ctx_boundary);
        read_neighbors(boundary, j_neighbors, n_neighbors, &neighbors);

        // the reciprocal volume calculation doesn't exactly fit typical SPH, so we'll implement it manually
        float recip_volume = 0;
        for(int k = 0; k < n_neighbors; k++){
            struct particle boundary_j = particle_at(&neighbors, k);
            
            recip_volume += W(euclid_dist(boundary[i].x, boundary[i].y, boundary_j.x, boundary_j.y) / H);
        }

        float m_i = boundary[i].rho / recip_volume;

        boundary[i].m = m_i;
    }
}

void calculate_density(struct particle *fluid, struct particle *boundary, struct neighbors_context *ctx_fluid, 
    struct neighbors_context *ctx_boundary)
{
    static float ones[MAX_POSSIBLE_NEIGHBORS] = {1, 1, 1, 1, 1, 1, 1, 1, 
        1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 
        1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
    int j_neighbors[MAX_POSSIBLE_NEIGHBORS], n_neighbors;
    struct particle_neighbors neighbors;

    #pragma omp for
    for(int i = 0; i < ctx_fluid->n_particles; i++){
        // fluid contribution to density of fluid
        n_neighbors = find_neighbors(j_neighbors, fluid, fluid, i, ctx_fluid);
        read_neighbors(fluid, j_neighbors, n_neighbors, &neighbors);
        float density_fluid_fluid = sph(ones, fluid[i], &neighbors, MASS);

        // boundary contribution to density of fluid
        n_neighbors = find_neighbors(j_neighbors, fluid, boundary, i, ctx_boundary);
        read_neighbors(boundary, j_neighbors, n_neighbors, &neighbors);
        float density_fluid_boundary = sph(ones, fluid[i], &neighbors, MASS);

        fluid[i].rho = density_fluid_fluid+density_fluid_boundary;
    }
}

// "slighly compressible" SPH, or weakly-compressible SPH (WCSPH) expresses pressure as an explicit function of density,
//  not an implicit one needing an iterative solver. Sometimes doesn't use the true speed of sound in a fluid, and
//  rather uses some number high enough that the density doesn't vary by too many percents. That's also how I did it
//  here because I saw that using the true speed caused instability. See Monaghan 1994 or Monaghan 2005.
void calculate_particle_pressure(struct particle *particles, int n_particles){
    #pragma omp for
    for(int i = 0; i < n_particles; i++){
        const float B = C*C*RHO_0/7;
        float rho_ratio = particles[i].rho/RHO_0;
        float pressure_i = B*(pow(rho_ratio, 7)-1);
        particles[i].p = (pressure_i > 0)? pressure_i : 0;
    }
}

void calculate_accelerations(float *du_dt_fluid, float *dv_dt_fluid, struct particle *fluid, 
    struct particle *boundary, struct neighbors_context *ctx_fluid, struct neighbors_context *ctx_boundary, 
    float gravity_x, float gravity_y)
{
    int j_neighbors[MAX_POSSIBLE_NEIGHBORS], n_neighbors;
    float temp_i[MAX_POSSIBLE_NEIGHBORS]; // holds a scalar quantity per neighbor to take the gradient of
    struct particle_neighbors neighbors;
    
    #pragma omp for
    for(int i = 0; i < ctx_fluid->n_particles; i++){
        // get the fluid neighbors of fluid[i]
        n_neighbors = find_neighbors(j_neighbors, fluid, fluid, i, ctx_fluid);
        read_neighbors(fluid, j_neighbors, n_neighbors, &neighbors);
        
        for(int k = 0; k < n_neighbors; k++){
            struct particle fluid_j = particle_at(&neighbors, k);

            // compute parts of the momentum-conserving pressure from fluid neighbors
            float pressure_ij = -(fluid[i].p/(fluid[i].rho*fluid[i].rho) + fluid_j.p/(fluid_j.rho*fluid_j.rho));
            
            // compute parts of the artificial pressure mentioned by Macklin 2013 (PBF) from fluid neighbors
            float W_ij = W(euclid_dist(fluid[i].x, fluid[i].y, fluid_j.x, fluid_j.y) / H);
            float artifical_pressure_ij = -0.1*powf(W_ij/W(0.2*H), 4);
            
            // compute parts of the viscosity from fluid neighbors
            float u_ij = fluid[i].u-fluid_j.u, v_ij = fluid[i].v-fluid_j.v;
            float x_ij = fluid[i].x-fluid_j.x, y_ij = fluid[i].y-fluid_j.y;
            float xy_dot_uv = x_ij*u_ij+y_ij*v_ij;
            float xy_dot_xy = x_ij*x_ij+y_ij*y_ij;
            float mu_ij = H*xy_dot_uv/(xy_dot_xy+0.01*H*H);
            float mean_rho = (fluid[i].rho+fluid_j.rho)/2;
            float viscosity_ij = 0.01*C*mu_ij/mean_rho;

            temp_i[k] = pressure_ij+artifical_pressure_ij+viscosity_ij;
        }

        // compute the acceleration due to fluid neighbors
        float2 acc_fluid_fluid_i = sph_gradient(temp_i, fluid[i], &neighbors, MASS);

        // get the boundary neighbors of fluid[i]
        n_neighbors = find_neighbors(j_neighbors, fluid, boundary, i, ctx_boundary);
        read_neighbors(boundary, j_neighbors, n_neighbors, &neighbors);
        
        for(int k = 0; k < n_neighbors; k++){
            struct particle boundary_j = particle_at(&neighbors, k);

            // compute parts of the momentum-conserving pressure from boundary neighbors
            float pressure_ij = -fluid[i].p/(fluid[i].rho*fluid[i].rho);
            
            // compute parts of the artificial pressure mentioned by Macklin 2013 (PBF) from boundary neighbors
            float W_ij = W(euclid_dist(fluid[i].x, fluid[i].y, boundary_j.x, boundary_j.y) / H);
            float artifical_pressure_ij = -0.1*powf(W_ij/W(0.2*H), 4);
            
            // compute parts of the viscosity from boundary neighbors
            float u_ij = fluid[i].u-boundary_j.u, v_ij = fluid[i].v-boundary_j.v;
            float x_ij = fluid[i].x-boundary_j.x, y_ij = fluid[i].y-boundary_j.y;
            float xy_dot_uv = x_ij*u_ij+y_ij*v_ij;
            float xy_dot_xy = x_ij*x_ij+y_ij*y_ij;
            float mu_ij = H*xy_dot_uv/(xy_dot_xy+0.01*H*H);
            float viscosity_ij = 0.01*C*mu_ij/fluid[i].rho; // use fluid density only, not the mean density

            temp_i[k] = pressure_ij+artifical_pressure_ij+viscosity_ij;
        }

        // compute the acceleration due to boundary neighbors
        float2 acc_fluid_boundary_i = sph_gradient(temp_i, fluid[i], &neighbors, MASS);

        du_dt_fluid[i] = gravity_x+acc_fluid_fluid_i.x+acc_fluid_boundary_i.x;
        dv_dt_fluid[i] = gravity_y+acc_fluid_fluid_i.y+acc_fluid_boundary_i.y;
    }
}


// METABALLS
// See Wikipedia article and original Blinn 1982 paper on metaballs. This function follows from the original derivation 
//  of a surface from points. Here I chose to hook off of the existing neighbor-finding code and kernel function.

void draw_metaballs(unsigned char *draw_buffer, struct particle *pixel_pseudoparticles, struct particle *fluid, 
    struct neighbors_context *ctx_fluid)
{
    int j_contributors[MAX_POSSIBLE_NEIGHBORS], n_contributors;
    struct particle_neighbors contributors;
    
    #pragma omp for collapse(2)
    for(int i = 0; i < 64; i++){
        for(int j = 0; j < 128; j++){
            int ij = i*128+j;

            n_contributors = find_neighbors(j_contributors, pixel_pseudoparticles, fluid, ij, ctx_fluid);
            read_neighbors(fluid, j_contributors, n_contributors, &contributors);

            float metaball_condition = 0;
            for(int k = 0; k < n_contributors; k++){
                struct particle fluid_j = particle_at(&contributors, k);

                float d = euclid_dist(pixel_pseudoparticles[ij].x, pixel_pseudoparticles[ij].y, fluid_j.x, fluid_j.y);

                // for now, using the SPH kernel (not normalized) as the Gaussian-like function suggested by Blinn
                float unnormalizing_factor = (4*M_PI*H*H)/7;
                metaball_condition += unnormalizing_factor*W(d/H);

                if(metaball_condition >= 1) break; // no need to keep adding contributions
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
void get_gravity(float2 *g){
    #ifdef MPU6050
    // We don't read in_accel_z_raw, which is the direction normal to the screen, but use in_accel_x_raw and 
    //  in_accel_y_raw as they are. This effectively implements a linear projection of the recorded gravity onto 
    //  the screen plane.
    int accel_x_raw = read_file_as_integer("/sys/bus/iio/devices/iio:device0/in_accel_x_raw"),
        accel_y_raw = read_file_as_integer("/sys/bus/iio/devices/iio:device0/in_accel_y_raw");
    // transform the raw values into a vector of magnitude being G at most
    g->x = (float)accel_y_raw / (1 << 14) * G; // (1 << 14) is about the vector magnitude reported when under gravity
    g->y = -(float)accel_x_raw / (1 << 14) * G;
    #else // if we don't have the MPU6050, just use a constant gravity vector
    g->x = 0;
    g->y = -G;
    #endif
}

void* get_gravity_routine(void *arg){
    float2 *g = (float2*)arg;
    
    struct timespec last, now;
    clock_gettime(CLOCK_MONOTONIC, &last);
    clock_gettime(CLOCK_MONOTONIC, &now);
    while(1){
        while(now.tv_nsec-last.tv_nsec < 1000000000/20){
            usleep(1000);
            clock_gettime(CLOCK_MONOTONIC, &now);
        }

        get_gravity(g);
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
    int particle_counter; // when the number of particles is determined numerically, this is used to count them


    int n_fluid;
    struct particle *fluid;
    float *du_dt, *dv_dt;
    
    // count the number of fluid particles we need
    n_fluid = 0;
    for(float x_0 = 0; x_0 < WIDTH; x_0 += R)
        for(float y_0 = 0; y_0 < HEIGHT; y_0 += R)
            if(in_initial_shape(x_0, y_0)) n_fluid++;
    
    // alloc fluid and derivatives
    fluid = (struct particle*)malloc(n_fluid*sizeof(struct particle));
    du_dt = (float*)malloc(n_fluid*sizeof(float));
    dv_dt = (float*)malloc(n_fluid*sizeof(float));
    
    // initialize fluid particles
    particle_counter = 0;
    for(float x_0 = 0; x_0 < WIDTH; x_0 += R){
        for(float y_0 = 0; y_0 < HEIGHT; y_0 += R){
            if(in_initial_shape(x_0, y_0)){
                fluid[particle_counter] = (struct particle){
                    .x = x_0, .y = y_0, .u = 0, .v = 0, 
                    .m = RHO_0*V, .rho = RHO_0 }; // no need to initialize .p (updated later)
                particle_counter++;
            }
        }
    }


    int n_boundary;
    struct particle *boundary;
    
    // count the number of boundary particles we need
    n_boundary = 0;
    for(float x_0 = 0; x_0 < WIDTH; x_0 += R) n_boundary += 2;
    for(float y_0 = 0; y_0 < HEIGHT; y_0 += R) n_boundary += 2;

    // alloc boundary particles
    boundary = (struct particle*)malloc(n_boundary*sizeof(struct particle));

    // initialize boundary particles properties
    particle_counter = 0;
    for(float x_0 = 0; x_0 < WIDTH; x_0 += R){
        boundary[particle_counter] = (struct particle){
            .x = x_0, .y = 0, .u = 0, .v = 0,
            .rho = RHO_0 }; // no need to initialize .m (calculated once later) or .p (never relevant)
        boundary[particle_counter+1] = (struct particle){
            .x = x_0, .y = HEIGHT, .u = 0, .v = 0,
            .rho = RHO_0 };
        particle_counter += 2;
    }
    for(float y_0 = 0; y_0 < HEIGHT; y_0 += R){
        boundary[particle_counter] = (struct particle){
            .x = 0, .y = y_0, .u = 0, .v = 0,
            .rho = RHO_0 };
        boundary[particle_counter+1] = (struct particle){
            .x = WIDTH, .y = y_0, .u = 0, .v = 0,
            .rho = RHO_0 };
        particle_counter += 2;
    }

    printf("dt = %f    (expected ticks/s) %d\n", DT, (int)(1/DT));
    printf("n_fluid = %d\n", n_fluid);
    printf("n_boundary = %d\n", n_boundary);


    struct timespec now; // initialize the time-keeping with the current time
    clock_gettime(CLOCK_MONOTONIC, &now);


    // initialize gravity and the gravity-reading thread
    float2 g;
    pthread_t gravity_thread;
    get_gravity(&g);
    pthread_create(&gravity_thread, NULL, get_gravity_routine, &g);


    struct particle *pixel_pseudoparticles;
    unsigned char *draw_buffer = (unsigned char*)calloc(1024, 1);
    struct timespec last_drew = now; // initialize the last-drew time to now
    pthread_t display_thread;

    // in leiu of defining a new function for finding the contributing particles to the metaballs condition, we'll just 
    //   reuse the neighbors search function (called with ctx_fluid as the argument)
    // to do so, we define pseudoparticles at the pixel centers (not unlike how we do treat the boundary)
    pixel_pseudoparticles = (struct particle*)malloc(64*128*sizeof(struct particle));
    for(int i = 0; i < 64; i++){
        for(int j = 0; j < 128; j++){
            float pixel_x = (j+0.5)*WIDTH/128, pixel_y = (64-(i+0.5))*HEIGHT/64;
            pixel_pseudoparticles[i*128+j].x = pixel_x;
            pixel_pseudoparticles[i*128+j].y = pixel_y;
        }
    }

    ssd1306_128x64_i2c_init(); // initialize display

    pthread_create(&display_thread, NULL, display_routine, draw_buffer); // launch the display thread


    // initialize statistics reporting
    float worst_max_rho_error_pct = 0, max_max_speed = 0;
    float t = 0, last_t = 0;
    struct timespec last_reported = now;


    // alloc neighbors search context
    const float x_min = 0-R, x_max = WIDTH+R, y_min = 0-R, y_max = HEIGHT+R;
    struct neighbors_context *ctx_fluid  = alloc_neighbors_context(n_fluid, x_min, x_max, y_min, y_max, 2*H), 
                             *ctx_boundary = alloc_neighbors_context(n_boundary, x_min, x_max, y_min, y_max, 2*H);

    // initialize context and masses of the boundary particles (these two functions never need to be called again)
    update_neighbors_context(ctx_boundary, boundary);
    calculate_boundary_pseudomass(boundary, ctx_boundary);

    // initialize du_dt and dv_dt with the accelerations for the zero-th time step
    update_neighbors_context(ctx_fluid, fluid);
    calculate_density(fluid, boundary, ctx_fluid, ctx_boundary);
    calculate_particle_pressure(fluid, n_fluid);
    calculate_accelerations(du_dt, dv_dt, fluid, boundary, ctx_fluid, ctx_boundary, g.x, g.y);


    // main loop, consisting of kick-drift-kick integration, drawing, and statistics reporting
    #pragma omp parallel num_threads(4)
    while(1){
        #pragma omp single
        {
            // kick: update the velocities by half a time step using the previous accelerations
            for(int i = 0; i < n_fluid; i++){
                fluid[i].u += 0.5*DT*du_dt[i];
                fluid[i].v += 0.5*DT*dv_dt[i];
            }

            // drift: update the positions by a full time step using the new velocities
            for(int i = 0; i < n_fluid; i++){
                fluid[i].x += DT*fluid[i].u;
                fluid[i].y += DT*fluid[i].v;
            }

            update_neighbors_context(ctx_fluid, fluid); // update the neighbors context for the new positions
        }

        // calculate the new accelerations from the new positions and new velocities
        calculate_density(fluid, boundary, ctx_fluid, ctx_boundary);
        calculate_particle_pressure(fluid, n_fluid);
        calculate_accelerations(du_dt, dv_dt, fluid, boundary, ctx_fluid, ctx_boundary, g.x, g.y);

        #pragma omp single
        {
            // kick: update the velocities by half a time step using the new accelerations
            for(int i = 0; i < n_fluid; i++){
                fluid[i].u += 0.5*DT*du_dt[i];
                fluid[i].v += 0.5*DT*dv_dt[i];
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
            float max_rho_error = 0, max_rho_error_pct;
            for(int i = 0; i < n_fluid; i++)
                if(fluid[i].rho > max_rho_error) max_rho_error = fluid[i].rho-RHO_0;
            max_rho_error_pct = max_rho_error/RHO_0*100;

            // compare avg_rho_ratio to worst_avg_rho_error_pct, which is the worst average rho ratio out of ALL frames
            // this is a critical statistic that shows to what degree the incompressibility constraint is being violated
            if(max_rho_error_pct > worst_max_rho_error_pct) worst_max_rho_error_pct = max_rho_error_pct;

            // record the maximum speed in a single frame as max_speed
            float max_speed = 0;
            for(int i = 0; i < n_fluid; i++){
                float speed = sqrt(fluid[i].u*fluid[i].u + fluid[i].v*fluid[i].v);
                if(speed > max_speed) max_speed = speed;
            }

            // compare max_speed to max_max_speed, which is the worst maximum speed out of ALL frames
            // this is a critical statistic that shows to what degree the maximum speed constraint is being violated
            if(max_speed > max_max_speed) max_max_speed = max_speed;
            
            // report frame rate and other statistics
            t += DT;
            if(t-last_t > 0.1){
                float elapsed = (now.tv_sec-last_reported.tv_sec) + (now.tv_nsec-last_reported.tv_nsec)/1e9;
                int tps = ((t-last_t)/DT)/elapsed;

                printf("sim time: %.2f, ", t);
                printf("ticks/s: %d, ", tps);
                printf("max rho error: %.3f\%% (worst) %.3f\%%, ", max_rho_error_pct, worst_max_rho_error_pct);
                printf("max speed: %.1f m/s (worst) %.1f m/s, ", max_speed, max_max_speed);
                printf("\n");

                last_t = t;
                last_reported = now;
            }
        }
    }
}