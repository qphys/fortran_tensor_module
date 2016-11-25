FC=gfortran
FFLAGS=-Ofast -fopenmp -Wall -fPIC
INCLUDES=./

all: utilities tensors
	
utilities:
	$(FC) -c -Wall ./utilities.f90
tensors:
	$(FC) -c $(FFLAGS)  ./tensors.f90
clean: 
	rm utilities.mod tensors.mod utilities.o tensors.o
	
	
