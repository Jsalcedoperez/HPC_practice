#include <stdio.h>
#include <mpi.h>
#include <stdlib.h>
int main(int argc, char *argv[])
{

  int ierr, procid, numprocs;
  
  ierr = MPI_Init(&argc, &argv);
  ierr = MPI_Comm_rank(MPI_COMM_WORLD, &procid);
  ierr = MPI_Comm_size(MPI_COMM_WORLD, &numprocs);


   double val = -1.0 * procid;
  //if (procid != 0){ 
   MPI_Send(&val, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
   printf("ProcID: %d\nsent value: %lf\nto Proc: 0\n",procid,val);
  //}
   if (procid == 0)
   {
   double val;
   double cum = 0.0;
   MPI_Status status;
   int i;
   for(i=0;i !=numprocs; ++i){
    
       ierr = MPI_Recv(&val,1,MPI_DOUBLE, MPI_ANY_SOURCE, 0 , MPI_COMM_WORLD,&status);
   
   if (ierr == MPI_SUCCESS){
      printf("ProcID: %d received value %lf.\n", procid,val);
      cum = cum + val;
   }
   else MPI_Abort(MPI_COMM_WORLD,1);
}

printf("the cumulative is: %d\n", cum);
}


ierr = MPI_Finalize();
}
