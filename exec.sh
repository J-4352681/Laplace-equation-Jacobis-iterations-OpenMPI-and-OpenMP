#!/bin/bash

arg4=$4
procs=${arg4:-1}

arg5=$5
export OMP_NUM_THREADS=${arg5:-1}

tpn=1
if [ "$1" == "mpi" ]; then
  tpn=8
fi

sbatch \
  -N $procs \
  --open-mode=append \
  --ntasks-per-node=$tpn \
  --partition=Blade \
  --exclusive \
  -o ./$1/exec/out.csv \
  -e ./$1/exec/err.txt \
  ./$1/run.sh ./$1/lplc-jcb-$1 $2 $3 $arg4 $arg5
