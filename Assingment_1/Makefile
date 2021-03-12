default:
	gfortran -fdefault-real-8 -fdefault-double-8 -o assignment_1 init.f90 physics.f90 derivs.f90 evolution.f90 assignment_1.f90
	./assignment_1

debug:
	gfortran -fdefault-real-8 -fdefault-double-8 -fcheck=all -Wall -o assignment_1 init.f90 physics.f90 derivs.f90 evolution.f90 assignment_1.f90
	./assignment_1

single : # single precision
	gfortran -fcheck=all -Wall -o assignment_1 init.f90 derivs.f90 assignment_1.f90
	./assignment_1
clean:
	rm *.out *.modd