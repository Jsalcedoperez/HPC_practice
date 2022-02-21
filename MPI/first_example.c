#include <stdio.h>
#include <mpi.h>

int main(int argc, char *argv[])
{

  int ierr, procid, numprocs;
  
  ierr = MPI_Init(&argc, &argv);
  ierr = MPI_Comm_rank(MPI_COMM_WORLD, &procid);
  ierr = MPI_Comm_size(MPI_COMM_WORLD, &numprocs);

  if (numprocs != 2)
  {

    printf("Error: number of procs is not 2!\n");
    return MPI_Abort(MPI_COMM_WORLD,1);
  }

  if (procid == 0)

  {

   //procid 0 will send number 3.14 to proc 1...
   double pi = 3.14;
   MPI_Send(&pi, 1, MPI_DOUBLE, 1, 0, MPI_COMM_WORLD);
   printf("ProcID: %d\nsent value: %lf\nto Proc: 1\n",procid,pi);}
  else
{
   //procid 1 will wait to receive a double from procid 0...
   double value;

   MPI_Status status;
   ierr = MPI_Recv(&value, 1,MPI_DOUBLE, 0, 0 , MPI_COMM_WORLD,&status);
   if (ierr == MPI_SUCCESS){
      printf("ProcID: %d received value %lf.\n", procid,value);
   }
   else{
      printf("ProcID %d did not received a value successfully!\n",procid);
}
}
    ierr = MPI_Finalize();


return 0;

 }
