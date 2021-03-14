default:
	gfortran -fdefault-real-8 -fdefault-double-8 -o main init.f90 outputs.f90 derivs.f90 evolution.f90 main.f90
	./main

iso:
	gfortran -fdefault-real-8 -fdefault-double-8 -o isothermal init.f90 outputs.f90 derivs.f90 evolution.f90 isothermal.f90
	./isothermal

sod:
	gfortran -fdefault-real-8 -fdefault-double-8 -o sod init.f90 outputs.f90 derivs.f90 evolution.f90 sod.f90
	./sod

debug:
	gfortran -fdefault-real-8 -fdefault-double-8 -fcheck=all -Wall -o main init.f90 outputs.f90 derivs.f90 evolution.f90 main.f90
	./main

single : # single precision
	gfortran -fcheck=all -Wall -o main init.f90 outputs.f90 derivs.f90 evolution.f90 main.f90
	./main
energy:
	gfortran -fdefault-real-8 -fdefault-double-8 -o main init.f90 outputs.f90 derivs.f90 evolution.f90 main.f90
	./main
	splash energy.out
clean:
	rm *.out *.mod output/snap*