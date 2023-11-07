/*
 * Parallel OpenMPI+OpenMP solution in C for bidimensional
 * Laplace equation using Jacobi's iterative method
*/

#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <math.h>
#include <mpi.h>
#include <string.h>
#include <omp.h>


#define COORD 0


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

int main(int argc, char *argv[]) {
  int procs, rank, provided;

  MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);
  MPI_Comm_size(MPI_COMM_WORLD, &procs);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  // Vars definition
  MPI_Status status1, status2;
  MPI_Request req1, req2;

  int i, r, c, rc, rn_up, rn_down;
  double *new, *old;
  int *senddispls, *sendcounts, *recvdispls, *recvcounts;
  double maxdiff = 0, diff, totalmaxdiff;



  if(argc != 3) {
    if(rank == 0) printf("Error en los parámetros. Uso: ./%s <matrix-size> <iterations>\n", argv[0]);
    exit(1);
  }

  const int mtrxsize = atoi(argv[1]);
  const int iters = floor(atoi(argv[2]) * 0.5);
  const int thds = atoi(getenv("OMP_NUM_THREADS"));

  if(mtrxsize <= 0 || iters <= 0 || procs <= 0 || thds <= 0) {
    if(rank == 0) printf("Los parámetros deben ser mayores a 0.\nEl tamaño de matriz debe ser divisible por la cantidad de hilos.\n", argv[0]);
    exit(1);
  }

  omp_set_num_threads(thds);



  const int lastproc = procs-1;

  const int lastcomm = iters*6+2;
  double commtime = 0, totaltime = 0, comms[lastcomm+3];
  const int timetaker = (procs > 1)?1:COORD;

  const int mtrxsize2 = mtrxsize*mtrxsize;

  const int useful = mtrxsize-2;
  const int mod = procs - (useful%procs);

  const int strip1 = useful/procs;
  const int strip2 = useful/procs+1;
  const int strip = (rank < mod)? strip1 : strip2;

  const int scattersize1 = strip1*mtrxsize;
  const int scattersize2 = strip2*mtrxsize;
  const int scattersize = (rank < mod)? scattersize1 : scattersize2;

  const int local1 = scattersize1 + 2*mtrxsize;
  const int local2 = scattersize2 + 2*mtrxsize;
  const int local = (rank < mod)? local1 : local2;

  const int localspace = local*sizeof(double);
  const int usefulmtrxsizespace = (mtrxsize-2)*sizeof(double);

  const int prev = rank-1, next = rank+1;

  const int colend = mtrxsize-1;

  const int lastrow = local - mtrxsize, penultrow = lastrow - mtrxsize;
  const int secondrow = mtrxsize, thirdrow = secondrow+mtrxsize;

  const int firstcolsecondrow = secondrow+1;
  const int lastcolsecondrow = secondrow+colend;
  const int firstcolpenultrow = penultrow+1;
  const int lastcolpenultrow = penultrow+colend;

  const int firstcollastrow = lastrow+1;


  // Init matrices
  if(rank == COORD) {
    senddispls = (int*)calloc(procs, sizeof(int));
    sendcounts = (int*)calloc(procs, sizeof(int));
    recvdispls = (int*)calloc(procs, sizeof(int));
    recvcounts = (int*)calloc(procs, sizeof(int));
    for (i = 0; i < procs; i++) {
      senddispls[i] = (i < mod)? scattersize1*i : scattersize1*mod + scattersize2*(i-mod);
      sendcounts[i] = (i < mod)? local1 : local2;
      recvcounts[i] = (i < mod)? scattersize1 : scattersize2;
      recvdispls[i] = senddispls[i] + mtrxsize;
    }

    new = (double*)calloc(mtrxsize2, sizeof(double));
    old = (double*)calloc(mtrxsize2, sizeof(double));
    init_matrix(old, mtrxsize);
    memcpy(new, old, mtrxsize2*sizeof(double));
  } else {
    new = (double*)calloc(local, sizeof(double));
    old = (double*)calloc(local, sizeof(double));
  }



  MPI_Barrier(MPI_COMM_WORLD); // Start processing


  if (rank==timetaker) comms[0] = MPI_Wtime();
  MPI_Scatterv(old, sendcounts, senddispls, MPI_DOUBLE, old, local, MPI_DOUBLE, COORD, MPI_COMM_WORLD);
  if (rank==timetaker) comms[1] = MPI_Wtime();

  if(rank != COORD) memcpy(new, old, localspace);


  #pragma omp parallel private(i,r,c,rc,rn_up,rn_down,diff)
  {


    for (i = 0; i < iters*6; i+=6) {

      #pragma omp for schedule(static)
      for (r = secondrow; r < lastrow; r += mtrxsize) {
        rn_up = r-mtrxsize;
        rn_down = r+mtrxsize;

        for (c = 1; c < colend; c ++) {
          rc = r+c;
          new[rc] = (
            old[rn_up+c] + old[rn_down+c] + old[rc-1] + old[rc+1]
          ) * 0.25;
        }
      }


      #pragma omp single nowait
      {
        if (rank==timetaker) comms[i+2] = MPI_Wtime();
        if(rank != 0) MPI_Isend(&new[secondrow],mtrxsize,MPI_DOUBLE,prev,2,MPI_COMM_WORLD,&req1);
        if(rank != lastproc) MPI_Isend(&new[penultrow],mtrxsize,MPI_DOUBLE,next,3,MPI_COMM_WORLD,&req2);
        if (rank==timetaker) comms[i+3] = MPI_Wtime();
      }


      #pragma omp for schedule(static)
      for (r = thirdrow; r < penultrow; r += mtrxsize) {
        rn_up = r-mtrxsize;
        rn_down = r+mtrxsize;

        for (c = 1; c < colend; c ++) {
          rc = r+c;
          old[rc] = (
            new[rn_up+c] + new[rn_down+c] + new[rc-1] + new[rc+1]
          ) * 0.25;
        }
      }


      #pragma omp single
      {
        if (rank==timetaker) comms[i+4] = MPI_Wtime();
        if(rank != lastproc) MPI_Recv(&new[lastrow],mtrxsize,MPI_DOUBLE,next,2,MPI_COMM_WORLD,&status1);
        if(rank != 0) MPI_Recv(new,mtrxsize,MPI_DOUBLE,prev,3,MPI_COMM_WORLD,&status2);
        if (rank==timetaker) comms[i+5] = MPI_Wtime();

        if(rank != 0) memcpy(&old[1], &new[1], usefulmtrxsizespace);
        if(rank != lastproc) memcpy(&old[firstcollastrow], &new[firstcollastrow], usefulmtrxsizespace);
      }


      #pragma omp for nowait schedule(static)
      for (c = firstcolsecondrow; c < lastcolsecondrow; c ++) {
        old[c] = (
          new[c-mtrxsize] + new[c+mtrxsize] + new[c-1] + new[c+1]
        ) * 0.25;
      }

      #pragma omp for nowait schedule(static)
      for (c = firstcolpenultrow; c < lastcolpenultrow; c ++) {
        old[c] = (
          new[c-mtrxsize] + new[c+mtrxsize] + new[c-1] + new[c+1]
        ) * 0.25;
      }



      // Asegura no modificar el buffer antes de que sea recibido
      #pragma omp single
      {
        if (rank==timetaker) comms[i+6] = MPI_Wtime();
        if(rank != lastproc) MPI_Wait(&req2,&status1);
        if(rank != 0) MPI_Wait(&req1,&status2);
        if (rank==timetaker) comms[i+7] = MPI_Wtime();
      }


    } // End iters



    #pragma omp for nowait reduction(max:maxdiff) schedule(static)
    for (r = secondrow; r < lastrow; r += mtrxsize) {
      for (c = 1; c < colend; c ++) {
        rc = r+c;
        diff = old[rc]-new[rc];
        if(maxdiff < diff) maxdiff = diff;
      }
    }


  } // End parallel



  if (rank==timetaker) comms[lastcomm] = MPI_Wtime();
  MPI_Gatherv(&old[secondrow], scattersize, MPI_DOUBLE, old, recvcounts, recvdispls, MPI_DOUBLE, COORD, MPI_COMM_WORLD);

  MPI_Reduce(&maxdiff, &totalmaxdiff, 1, MPI_DOUBLE, MPI_MAX, COORD, MPI_COMM_WORLD);
  if (rank==timetaker) comms[lastcomm+1] = MPI_Wtime(); // End processing



  if(rank==timetaker && procs > 1) MPI_Send(comms,lastcomm+2,MPI_DOUBLE,COORD,0,MPI_COMM_WORLD);
  if(rank==COORD && procs > 1) MPI_Recv(comms,lastcomm+2,MPI_DOUBLE,timetaker,0,MPI_COMM_WORLD,&status1);


  if(rank == COORD) {
    totaltime = comms[lastcomm+1] - comms[0];

    for(r = 0; r < lastcomm+2; r+=2) {
      commtime += comms[r+1] - comms[r];
    }

    // print_float_matrix(old, mtrxsize);
    // print_float_matrix(new, mtrxsize);

    printf(
      "%i,%i,%i,%i,%f,%f,%f\n",
      mtrxsize, iters*2, procs, thds, totalmaxdiff, totaltime, commtime
    );
  }

  free(new); free(old);

  MPI_Finalize();

  return 0;
}
