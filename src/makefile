COMPILER= gfortran

FLAGS = -O3 #-fopenmp

EXEC = poisson

SRC = $(wildcard *.f90) 

OBJ = $(SRC:.f90=.o)

$(EXEC): $(OBJ)
	$(COMPILER) $(FLAGS) -o $@ $^ $(LAPACK)

types.o: types.f90
	$(COMPILER) $(FLAGS) -c $<

quadrature.o: quadrature.f90 types.o
	$(COMPILER) $(FLAGS) -c $<

read_write.o: read_write.f90 types.o potential.o
	$(COMPILER) $(FLAGS) -c $<

potential.o: potential.f90 types.o quadrature.o
	$(COMPILER) $(FLAGS) -c $<

main.o: main.f90 types.o potential.o 	read_write.o 
	$(COMPILER) $(FLAGS) -c $<

clean:
	rm -rf *.o *.mod

mrproper: clean
	rm -rf $(EXEC)
