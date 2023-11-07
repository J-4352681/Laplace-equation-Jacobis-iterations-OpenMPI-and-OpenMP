## Bidimensional Laplace Ecuation with Jacobi iterative method

### Sequential solution written in C

- Dependencies:
  - C programming langugage compiler and executor
  - OpenMP for C
  - OpenMPI for C
  - Slurm (for task scheduling)

- Scripts execution:

  1. Open a console in this directory.
  2. Compile the program of every subdirectory with `./compile.sh`.
  3. Call `./exec.sh` script the following way:
    
    `./exec.sh <mpi | mpiomp | seq> <lado de matriz> <iteraciones> [<nodos> [<hilos>]]`
    
    - If it's executed with "mpi", the amount of processors will be 8 by node.
    - If it's exectuted witn "mpiomp", the amount of processors will be equal to the amount of nodes passed in the arguments.
    - If the amount of threads isn't assigne, it'll be equal to 1 (useless with "mpi").
  4. The results will be saved in the "out.csv" file of the corresponding directory.
