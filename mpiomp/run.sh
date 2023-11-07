#!/bin/bash
export OMP_NUM_THREADS=$5
mpirun --bind-to none $1 $2 $3 
