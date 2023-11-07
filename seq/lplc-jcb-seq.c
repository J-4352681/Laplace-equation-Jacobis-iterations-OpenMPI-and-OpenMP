#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <math.h>
#include <string.h>
#include <float.h>

double dwalltime(){
  struct timeval tv;
  gettimeofday(&tv,NULL);
  return tv.tv_sec + tv.tv_usec/1000000.0;
}

void print_float_matrix(double *mtrx, int size) {
  int i, j, in;
  printf("\n");
  for(i=0; i<size; i++) {
    in=i*size;
    printf("| ");
    for(j=0; j<size; j++) {
      printf("%.3f | ", mtrx[in+j]);
    }
    printf("\n");
  }
}

void init_matrix(double *mtrx, int mtrxsize) {
  int i;
  double init_values[] = {3.5, 2.4, 7.5, 9.4, 2.5, 4.2};
  // double init_values[] = {2.7};
  int size_values = sizeof(init_values) / sizeof(init_values[0]);
  for(i = 0; i < mtrxsize; i++)
    mtrx[i] = init_values[i%size_values];
  for(i = 1; i < mtrxsize-1; i++)
    mtrx[i*mtrxsize] = init_values[i%size_values];
  for(i = 1; i < mtrxsize-1; i++)
    mtrx[i*mtrxsize+(mtrxsize-1)] = init_values[i%size_values];
  for(i = 0; i < mtrxsize; i++)
    mtrx[i+(mtrxsize*(mtrxsize-1))] = init_values[i%size_values];
}

int main(int argc, char *argv[]) {
  int mtrxsize, mtrxsize2, iters, rowend, colend, rowstart;

  if(argc != 3) {
    printf("Error en los parámetros. Uso: ./%s <matrix-size> <iterations>\n", argv[0]);
    exit(1);
  }

  mtrxsize = atoi(argv[1]);
  iters = atoi(argv[2]);

  colend = mtrxsize-2;
  mtrxsize2 = mtrxsize*mtrxsize;
  rowend = mtrxsize2-mtrxsize;
  rowstart = mtrxsize+1;

  int i, r, c, rc, rn_up = 1, rn_down;

  double time;

  double *new, *old;
  new = (double*)calloc(mtrxsize2, sizeof(double));
  old = (double*)calloc(mtrxsize2, sizeof(double));
  double maxdiff = 0, diff;
  init_matrix(old, mtrxsize);
  init_matrix(new, mtrxsize);


  time = dwalltime();

  for (i = 0; i < iters*0.5; i++) {
    for (r = rowstart; r < rowend; r += mtrxsize) {
      rn_up = r-mtrxsize;
      rn_down = r+mtrxsize;

      for (c = 0; c < colend; c ++) {
        rc = r+c;
        new[rc] = (
          old[rn_up+c] + old[rn_down+c] + old[rc-1] + old[rc+1]
        ) * 0.25;
      }
    }

    for (r = rowstart; r < rowend; r += mtrxsize) {
      rn_up = r-mtrxsize;
      rn_down = r+mtrxsize;

      for (c = 0; c < colend; c ++) {
        rc = r+c;
        old[rc] = (
          new[rn_up+c] + new[rn_down+c] + new[rc-1] + new[rc+1]
        ) * 0.25;
      }
    }
  }

  for (r = rowstart; r < rowend; r += mtrxsize) {
    for (c = 0; c < colend; c ++) {
      rc = r+c;
      diff = old[rc]-new[rc];
      if(maxdiff < diff) maxdiff = diff;
    }
  }

  time = dwalltime()-time;
  
  
  
  // print_float_matrix(new, mtrxsize);
  // print_float_matrix(old, mtrxsize);
  
  printf("%i,%i,%f,%f\n", mtrxsize, iters, maxdiff, time);

  return 0;
}
