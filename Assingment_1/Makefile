FFLAGS= -fdefault-real-8 -fdefault-double-8 -O3
FILES= init.f90 outputs.f90 derivs.f90 evolution.f90

default:
	# change defaults
	cp wave.splash.defaults splash.defaults
	gfortran ${FFLAGS} -o main ${FILES} main.f90
	./main

iso:
	# change defaults
	cp iso.splash.defaults splash.defaults
	gfortran ${FFLAGS} -o isothermal ${FILES} isothermal.f90
	./isothermal

sod:
	# change defaults
	cp sod.splash.defaults splash.defaults 
	gfortran ${FFLAGS} -o sod ${FILES} sod.f90
	./sod

debug:
	gfortran ${FFLAGS} -fcheck=all -Wall -o main ${FILES} main.f90
	./main

single: # single precision
	gfortran -fcheck=all -Wall -o main ${FILES} main.f90
	./main
energy:
	cp energy.splash.defaults splash.defaults 
	gfortran ${FFLAGS} -o main ${FILES} main.f90
	./main
clean:
	rm *.out *.mod output/snap* main isothermal sod