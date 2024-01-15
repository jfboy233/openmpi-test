#!/bin/bash

#SBATCH -N 1
#SBATCH -n 6

mpiexec --allow-run-as-root ./helloworld
# echo helloworld