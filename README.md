# pi-sph-fluid

![headliner](headliner.jpg)

This project uses an SSD1306 OLED module, an MPU6050 accelerometer, and a Raspberry Pi 4 to realize a fluid simulation toy that lets you toss around an ocean of water in your hand. (That is, it simulates fluid dynamics at a scale larger than the actual screen.) You can see it in action on [Youtube](https://youtu.be/6_YntT-Zha0). This can be thought of as a "free-surface flow" problem, and it's solved by a "smoothed particle hydrodynamics" (SPH) technique.

## Running

The only dependencies are OpenMP and the [ssd1306](https://github.com/lexus2k/ssd1306) driver by lexus2k. Thanks to the SDL emulation of the SSD1306 offered by the ssd1306 library, this project offers a choice of two executables: `desktop_sph_fluid` and `pi_sph_fluid`.

1. `desktop_sph_fluid` generates and renders a circular drop falling on a dry surface, and it can be compiled (on Linux) with the commands

```bash
git clone https://github.com/colonelwatch/pi-sph-fluid
cd pi-sph-fluid
git submodule update --init ssd1306
make desktop_sph_fluid

./desktop_sph_fluid
```

2. `pi_sph_fluid` generates the same but with gravity determined by reading the MPU6050 accelerometer. I2C and the driver for accelerometer need to be enabled by configuring `/boot/config.txt`. All the necessary commands are

```bash
sudo apt install git

# configure /boot/config.txt (you can also write these lines in manually)
sudo raspi-config nonint do_i2c 0 # call raspi-config non-interactively for turning on i2c correctly
sudo raspi-config nonint set_config_var dtparam=i2c_arm on,1000000 /boot/config.txt # set i2c speed to 1MHz
echo dtoverlay=mpu6050 | sudo tee -a /boot/config.txt # pull up device tree overlay for mpu6050
sudo reboot

git clone https://github.com/colonelwatch/pi-sph-fluid
cd pi-sph-fluid
git submodule update --init ssd1306
make pi_sph_fluid

# limit the rate at which the MPU6050 is being polled
echo 10 | sudo tee /sys/bus/iio/devices/iio:device0/sampling_frequency

./pi_sph_fluid
```

## What's implemented?

Besides the ssd1306 driver, this project is just under 750 lines of C! Here's what this project implements in that many lines:

1. Linked-list neighbors search, a common technique but outlined well in [1]
2. Weakly-compressible smoothed particle hydrodynamics (WCSPH), completely described in [2] and [3]
    * Artificial viscosity and momentum-preserving pressure, described in the same
    * Ordinary advection (not XSPH)
    * Working around negative density error by clamping negative pressure to zero, mentioned in the IISPH, DFSPH, and PBF papers ([4-6] respectively) along with other papers
3. Boundary handling offered by [7]
4. Fake surface tension effects using artificial pressure, mentioned in [6]
5. Rendering using metaballs, following the original implementation in [8]
6. OpenMP acceleration

## What's not implemented?

This project really sits on the ground floor of SPH, and some next steps include:

1. Z-order or some other -order sort of particles for a higher cache hit rate, as described in [1]
    * Ideally, this sort should be multithreaded
2. Inferring velocity from accelerations taken from the MPU6050 and assigning that velocity to the boundary
    * This would add a bit more realism to the fluid simulation
3. GPU acceleration or an implementation of one of the incompressible schemes (PCISPH, IISPH, DFSPH, PBF, etc)
    * There's a failed attempt at IISPH in the `IISPH` branch
    * Some of the most interesting fluid phenomena don't arise in the small number of particles that CPU WCSPH can handle without going unstable (currently 650)

## References

[1] J. M. Domínguez, A. J. C. Crespo, M. Gómez-Gesteira, and J. C. Marongiu, “Neighbour lists in smoothed particle hydrodynamics,” International Journal for Numerical Methods in Fluids, vol. 67, no. 12, pp. 2026–2042, 2011, doi: 10.1002/fld.2481.

[2] J. J. Monaghan, “Smoothed particle hydrodynamics,” Reports on progress in physics, vol. 68, no. 8, p. 1703, 2005.

[3] J. J. Monaghan, “Simulating Free Surface Flows with SPH,” Journal of Computational Physics, vol. 110, no. 2, pp. 399–406, 1994, doi: 10.1006/jcph.1994.1034.

[4] M. Ihmsen, J. Cornelis, B. Solenthaler, C. Horvath, and M. Teschner, “Implicit Incompressible SPH,” IEEE Transactions on Visualization and Computer Graphics, vol. 20, no. 3, pp. 426–435, Mar. 2014, doi: 10.1109/TVCG.2013.105.

[5] J. Bender and D. Koschier, “Divergence-free smoothed particle hydrodynamics,” in Proceedings of the 14th ACM SIGGRAPH / Eurographics Symposium on Computer Animation, Los Angeles California: ACM, Aug. 2015, pp. 147–155. doi: 10.1145/2786784.2786796.

[6] M. Macklin and M. Müller, “Position based fluids,” ACM Trans. Graph., vol. 32, no. 4, pp. 1–12, Jul. 2013, doi: 10.1145/2461912.2461984.

[7] N. Akinci, M. Ihmsen, G. Akinci, B. Solenthaler, and M. Teschner, “Versatile rigid-fluid coupling for incompressible SPH,” ACM Trans. Graph., vol. 31, no. 4, pp. 1–8, Aug. 2012, doi: 10.1145/2185520.2185558.

[8] J. F. Blinn, “A Generalization of Algebraic Surface Drawing,” ACM Trans. Graph., vol. 1, no. 3, pp. 235–256, Jul. 1982, doi: 10.1145/357306.357310.