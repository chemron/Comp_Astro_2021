default:
	gfortran -fdefault-real-8 -fdefault-double-8 -fcheck=all -Wall -o assignment_1 init.f90 derivs.f90 assignment_1.f90
	./assignment_1
single : # single precision
	gfortran -fcheck=all -Wall -o assignment_1 init.f90 derivs.f90 assignment_1.f90
	./assignment_1
clean:
	rm *.out *.mod