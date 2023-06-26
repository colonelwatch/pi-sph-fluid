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
# configure /boot/config.txt (you can also write these lines in manually)
echo dtparam=i2c_arm=on,1000000 | sudo tee -a /boot/config.txt # enable i2c with a default speed of 1 MHz
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

1. Linked-list neighbors search, a common technique but outlined well in Domínguez 2011
2. Weakly-compressible smoothed particle hydrodynamics (WCSPH), completely described in Monaghan 2005 and Monaghan 1994
    * Artificial viscosity and momentum-preserving pressure, described in the same
    * Ordinary advection (not XSPH)
    * Working around negative density error by clamping negative pressure to zero, mentioned in the IISPH, DFSPH, and PBF papers (Ihmsen 2014, Bender 2015, and Macklin 2013 respectively) along with other papers
3. Boundary handling offered by Akinci 2012
4. Fake surface tension effects using artificial pressure, mentioned in Macklin 2013
5. Rendering using metaballs, following the original implementation in Blinn 1982
6. OpenMP acceleration

## What's not implemented?

This project really sits on the ground floor of SPH, and some next steps include:

1. Z-order or some other -order sort of particles for a higher cache hit rate
    * Ideally, this sort should be multithreaded
2. Inferring velocity from accelerations taken from the MPU6050 and assigning that velocity to the boundary
    * This would add a bit more realism to the fluid simulation
3. GPU acceleration or an implementation of one of the incompressible schemes (PCISPH, IISPH, DFSPH, PBF, etc)
    * There's a failed attempt at IISPH in the `IISPH` branch
    * Some of the most interesting fluid phenomena don't arise in the small number of particles that CPU WCSPH can handle without going unstable (currently 650)