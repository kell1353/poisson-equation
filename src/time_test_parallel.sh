#!/bin/bash
# shell script for running program with thread numbers

# Clear old outputs and compile the program
make mrproper
make

# Iterate over the number of threads
# and run the program for each thread
# number
for i in {1..4}
do
	# set number of threads to i
        export OMP_NUM_THREADS=$i   
        for j in {1..10}
        do
        # run executable with namelist for given thread number i	
        ./poisson box_parameters.namelist
    done
done

