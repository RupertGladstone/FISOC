#!/bin/bash
#PBS -q normal
#PBS -P m68
#PBS -lncpus=32
#PBS -lmem=3000mb
#PBS -lwalltime=00:20:00
#PBS -l wd

module load openmpi

mpirun -np 32 ElmerSolver_mpi 
