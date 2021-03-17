FFLAGS= -fdefault-real-8 -fdefault-double-8 -Ofast
FILES= init.f90 outputs.f90 derivs.f90 evolution.f90

default:
	gfortran ${FFLAGS} -o main ${FILES} main.f90
	./main

iso:
	gfortran ${FFLAGS} -o isothermal ${FILES} isothermal.f90
	./isothermal

sod:
	gfortran ${FFLAGS} -o sod ${FILES} sod.f90
	./sod

debug:
	gfortran ${FFLAGS} -fcheck=all -Wall -o main ${FILES} main.f90
	./main

single : # single precision
	gfortran -fcheck=all -Wall -o main ${FILES} main.f90
	./main
energy:
	gfortran ${FFLAGS} -o main ${FILES} main.f90
	./main
	splash energy.out
clean:
	rm *.out *.mod output/snap* main isothermal sod