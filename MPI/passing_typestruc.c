#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

void Build_derived_type(float* a_ptr, float* b_ptr, int* n_ptr, MPI_Datatype* mesg_mpi_t_ptr)
     {
     int block_lengths[3]; //each block only has one element because they are not contigious.
     MPI_Aint displacements[3];
     MPI_Datatype typelist[3];

     MPI_Aint start_address;
     MPI_Aint address;
     block_lengths[0] = block_lengths[1] = block_lengths[2] = 1;
     
     typelist[0] = MPI_FLOAT;
     typelist[1] = MPI_FLOAT;
     typelist[2] = MPI_INT;
     /*first element, a, is at displacement 0*/
     displacements[0] = 0;
     MPI_Address(a_ptr, &start_address);
     MPI_Address(b_ptr, &address);
     displacements[1] = address - start_address;

     MPI_Address(n_ptr, &address);
     displacements[2] = address - start_address;
     

     MPI_Type_struct(3,block_lengths, displacements, typelist, mesg_mpi_t_ptr);
     MPI_Type_commit(mesg_mpi_t_ptr);

}


void Get_data(float* a_ptr, float* b_ptr, int* n_ptr, int my_rank){
              if (my_rank == 0){
                printf("%f, %f, %d", a_ptr, b_ptr, n_ptr);
              }

     Build_derived_type(a_ptr,b_ptr,n_ptr,&mesg_mpi_t);
     MPI_Bcast(a_ptr,1,mesg_mpi_t,0,MPI_COMM_WORLD); //passing only pointer of first element.
     }

