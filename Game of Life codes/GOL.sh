#!/bin/bash
#
## SLURM submission script for MPI job  Game of Life
## change ntasks to change number of processes
#SBATCH --job-name=GOL
#SBATCH --output=GOL%J.out
#SBATCH --error=GOL%J.err
#
#SBATCH --ntasks=1
#SBATCH --time=5:00
#SBATCH --partition=class

## change n to modify number cells
## change g to modify number of generations
## set -debug 1 to see output of final population (only small n!)
mpiexec ./game -n 10 -g 5 -debug 0
