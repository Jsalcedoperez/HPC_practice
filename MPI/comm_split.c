#include <stdio.h>
#include <math.h>
#include "mpi.h"

int main(int argc, char* argv[]){

int p,my_rank;
MPI_Comm my_row_comm;
int my_row, my_rank_in_row;
int q, test;

MPI_Init(&argc,&argv);
MPI_Comm_size(MPI_COMM_WORLD,&p);
MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

q = (int) sqrt((double) p);
/* my_rank is rank in MPI_COMM_WORLD
   q * q = p */
my_row = my_rank / q;

MPI_Comm_split(MPI_COMM_WORLD,my_row,my_rank,&my_row_comm);

/* Test the new communicator */

MPI_Comm_rank(my_row_comm, &my_rank_in_row);
if (my_rank_in_row == 0){
  test = my_row;
}else{
   test = 0;
 }

//MPI_Bcast(&test, 1, MPI_INT, 0, my_row_comm);

printf("Process %d > my_row = %d, my_rank_in_row = %d, test = %d\n",
        my_rank, my_row,my_rank_in_row,test);
MPI_Finalize();
return 0;
}
/*
 *
 * MPI_Cart_sub(MPI_Comm comm, int* free_coords, MPI Comm *newcomm)
 * partitions a communicator into subgroups which form lower-dimensional Cartesian subgrids.
 * comm communicator with cartesian structure,
 * free_coords - an array which specifies which dimensions are free (true)
 * and which are not free (false)
 * newcomm - communicator containing the subgrid that includes the calling processs.
 *
 int free_coords[2];
 MPI_Comm_row_comm;
 
 free_coords[0] = 0;
 free_coords[1] = 1;
 MPI_Cart_sub(grid_comm, free_coords, &row_comm)


}
