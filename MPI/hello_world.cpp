#include <stdio.h>
#include "mpi.h"
#include <iostream>
int main(int argc, char *argv[])
{

  int ierr, procid, numprocs;
  ierr = MPI_Init(&argc,&argv);
  ierr = MPI_Comm_rank(MPI_COMM_WORLD,&procid);
  ierr = MPI_Comm_size(MPI_COMM_WORLD,&numprocs);

  printf("hello world! I am process %d out of %d!\n", procid,numprocs);
  /*cout << "hello world! I am process "<<procid<<"out of "<<numprocd<<endl;*/
  ierr = MPI_Finalize();
  return 0;

}
