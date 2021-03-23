FFLAGS= -fdefault-real-8 -fdefault-double-8 -O3
FILES= integrators.f90 output.f90 main.f90

default:
	mkdir -p output
	gfortran ${FFLAGS} -o main ${FILES}
	./main

debug:
	gfortran ${FFLAGS} -fcheck=all -Wall -fbacktrace -Wextra -o main ${FILES}

clean:
	rm *.out *.mod output/snap* main