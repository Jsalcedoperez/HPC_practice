#include <stdio.h>
#include <mpi.h>

int main(int argc, char *argv[])
{

  int ierr, procid, numprocs;
  
  ierr = MPI_Init(&argc, &argv);
  ierr = MPI_Comm_rank(MPI_COMM_WORLD, &procid);
  ierr = MPI_Comm_size(MPI_COMM_WORLD, &numprocs);

  if (procid%2 == 0)

  {

   //procid 0 will send number 3.14 to proc 1...
   double value = 3.14+procid;
   MPI_Send(&value, 1, MPI_DOUBLE, procid+1, 0, MPI_COMM_WORLD);
   printf("ProcID: %d\nsent value: %lf\nto Proc: %d\n",procid,value,procid+1);}
  else
{
   //procid 1 will wait to receive a double from procid 0...
   double val;

   MPI_Status status;
   ierr = MPI_Recv(&val,1,MPI_DOUBLE, MPI_ANY_SOURCE, 0 , MPI_COMM_WORLD,&status);
   if (ierr == MPI_SUCCESS){
      printf("ProcID: %d received value %lf.\n", procid,val);
   }
   else{
      printf("ProcID %d did not successfully receive a value!\n",procid);
}
}
    ierr = MPI_Finalize();

    return 0;

}
