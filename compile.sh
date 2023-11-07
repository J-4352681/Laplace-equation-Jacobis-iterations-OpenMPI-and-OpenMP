#!/bin/bash

gcc -O3 ./seq/lplc-jcb-seq.c -o ./seq/lplc-jcb-seq

mpicc -O3 ./mpi/lplc-jcb-mpi.c -o ./mpi/lplc-jcb-mpi -lm

mpicc -fopenmp -O3 ./mpiomp/lplc-jcb-mpiomp.c -o ./mpiomp/lplc-jcb-mpiomp -lm
