CC := gcc
CFLAGS := -Wall -Ofast -march=native
INCLUDES := -Iexternal/ssd1306/src
LDFLAGS := -lm -pthread -fopenmp
LIBS :=

.PHONY: all clean

all: pi_sph_fluid

clean:
	cd external/ssd1306/src && make -f Makefile.linux clean
	cd external/ssd1306/tools/sdl && make -f Makefile.linux clean
	rm -f pi_sph_fluid
	rm -f desktop_sph_fluid


desktop_sph_fluid: src/pi_sph_fluid.c
	cd external/ssd1306/src && make -f Makefile.linux ../bld/libssd1306.a SDL_EMULATION=y
	cd external/ssd1306/tools/sdl && make -f Makefile.linux ../../bld/libssd1306_sdl.a
	$(CC) $^ external/ssd1306/bld/libssd1306.a external/ssd1306/bld/libssd1306_sdl.a -o $@ \
		$(INCLUDES) -Iexternal/ssd1306/tools/sdl -I/usr/include/SDL2 $(CFLAGS) -DSDL_EMULATION \
		$(LDFLAGS) -lSDL2 $(LIBS) -L/usr/lib/libSDL2-2.0.so

pi_sph_fluid: src/pi_sph_fluid.c
	cd external/ssd1306/src && make -f Makefile.linux ../bld/libssd1306.a
	$(CC) $^ external/ssd1306/bld/libssd1306.a -o $@ $(INCLUDES) $(CFLAGS) -DMPU6050 $(LDFLAGS) $(LIBS)