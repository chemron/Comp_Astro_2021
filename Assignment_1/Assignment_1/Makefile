FFLAGS= -fdefault-real-8 -fdefault-double-8 -O3
FILES= init.f90 edges.f90 outputs.f90 derivs.f90 evolution.f90
DEBUG=sod

default:
	gfortran ${FFLAGS} -o main ${FILES} main.f90
	./main

iso:
	gfortran ${FFLAGS} -o isothermal ${FILES} isothermal.f90
	./isothermal

sod:
	gfortran ${FFLAGS} -o sod_shock ${FILES} sod_shock.f90
	./sod_shock