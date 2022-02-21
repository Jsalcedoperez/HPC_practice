#include <stdio.h>
#include "mpi.h"
#include <math.h>

int main(int argc, char* argv[]){

  int p, my_rank, row_test,q;
  MPI_Comm grid_comm;
  int dim_sizes[2];
  int wrap_around[2]; //indicates whether grid is periodic.
  int coordinates[2];
  int free_coords[2];
  int reorder = 1; //whether or not ranking of initial process may be reordered. 
  int my_grid_rank, grid_rank;
  MPI_Comm row_comm; //we are going to make subsets of the grid_communicator that go along the rows of the cartensian grid_comunicator.  
  MPI_Comm col_comm; //we are going to make subsets of the grid_communicator that go along the columns of the cartensian grid_comunicator. 
  MPI_Init(&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD, &p);
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

  q = (int) sqrt((double) p);
  dim_sizes[0] = dim_sizes[1] = q;
  wrap_around[0] = wrap_around[1] = 1;
  MPI_Cart_create(MPI_COMM_WORLD,2,dim_sizes,wrap_around,reorder,&grid_comm);
  MPI_Comm_rank(grid_comm,&my_grid_rank);
  MPI_Cart_coords(grid_comm,my_grid_rank,2,coordinates);
  //going from coordinates to rank.
  MPI_Cart_rank(grid_comm, coordinates, &grid_rank);
  printf("Process %d > my_grid_rank = %d,coords = (%d,%d), grid_rank = %d\n",my_rank,my_grid_rank,coordinates[0],coordinates[1],grid_rank);
  free_coords[0] = 0;
  free_coords[1] = 1;
  //partion a communicator into subgroups which form lower-dimensional cartensian subgrids. 
  MPI_Cart_sub(grid_comm,free_coords,&row_comm);

  if(coordinates[1] == 0){
   row_test = coordinates[0];}
  else{
   row_test = -1;} 

  //MPI_Bcast(&row_test, 1, MPI_INT, 0, row_comm);

  printf("Process %d > coords = (%d,%d), row_test = %d\n",my_rank,coordinates[0],coordinates[1], row_test);
  
  /*
    int MPI_Cart_shift(MPI_Comm comm, int direction, int displ, int *source, int *dest)
      comm-communicator with cartensian structure.
      direction - coordinate dimension of shift, in range [0,n-1] for an n-dimensional cartensian grid
      displ - displacement (>0: upwards shift, <0: downwards shift), with periodic wraparound possible if communicator created with periodic boundary conditions turned on.
      outputs are possible inputs to MPI_Sendrecv
      source -rank of process to receive data from, obtained by substracting displ from coordinate determined by direction.
      dest- rank of process to send data to, obtained by adding displ to coordinate determined direction.
      these may be undefined (i.e., = MPI_PROC_NULL) if shift points outside grid structure and the periodic boundary conditions off. 
      examples:
      MPI_Cart_shift(comm, 1,1, &source,&dest) --> moving along columns in forward direction
      MPI_Cart_shift(comm, 0,-1, &source,&dest) --> moving along rows in backward direction.
     
     
  */
      
  
  MPI_Finalize();return 0;}
