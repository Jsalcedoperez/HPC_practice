#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

/* func.c */

float f(float);
float Trap(float, float, int, float);
void Get_data(float*, float*, int*,
              int);
int main(int argc, char** argv)
{
int my_rank;        //procs ranks
int p;              //number of process
float a = 1.0;      //left end point
float b = 6.0;    //right end point
int n = 1024;         //number of traps;
float local_a ;      //left endpoint my process
float local_b;       //right endpoint my process
int local_n;          //number of trapezoids
float integral;       //integral over my interval
float total = -1;
int source;
int dest  = 0;
int tag = 0;
float h;

MPI_Status status;
MPI_Init(&argc,&argv);
MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
MPI_Comm_size(MPI_COMM_WORLD, &p);

h = (b-a) / n;  // trap based length
local_n = n/p;

local_a = a + my_rank*local_n*h;
local_b = local_a + local_n*h;

integral = Trap(local_a,local_b,local_n,h);

/*MPI_Reduce(void *sendbuf, void *recvbuf, int count,
            MPI Datatype datatype, MPI_Op op, int root,
            MPI_Comm comm);

Must be called in all processes in a communicator, but result only
available in root process
*/

MPI_Reduce(&integral, &total, 1, MPI_FLOAT,MPI_SUM,0,MPI_COMM_WORLD);
/* this is without using the reducing intrinsics.
if (my_rank == 0){
total = integral;
for (source = 1; source <p; source++){
  MPI_Recv(&integral, 1, MPI_FLOAT, source, tag, MPI_COMM_WORLD,&status);
  printf("PE %d <- %d, %f\n", my_rank,source,integral);
  total = total + integral;



 }
}

else{

   MPI_Send(&integral, 1, MPI_FLOAT, dest, tag, MPI_COMM_WORLD);
}

*/
if (my_rank == 0){

printf("with n = %d trapezoids, our estimate\n",n);
printf("of the integral from %f to %f = %f\n",a,b,total);
}
MPI_Finalize();
return 0;

}


float f(float x)
{

return x*x;

}


float Trap(float a, float b, int n, float h){

float integral;
float x;
int i;

integral = (f(a) + f(b))/2.0;
x = a;
for (i=1; i<= n-1; i++)
{
x = x + h;
integral = integral + f(x);
}
return integral*h;
}
void Get_data(float* a_ptr, float* b_ptr, int* n_ptr,
              int my_rank){

      if (my_rank == 0){
          printf("Enter a, b and n\n");
          scanf("%f %f %d", a_ptr, b_ptr, n_ptr);
      }
      /*
        int MPI_Bcast(void *buffer, int count, MPI_Datatype datatype, int root, MPI_Comm comm)
      MPI_Bcast(a_ptr,1,MPI_FLOAT,0,MPI_COMM_WORLD);
      MPI_Bcast(b_ptr,1,MPI_FLOAT,0,MPI_COMM_WORLD);
      MPI_Bcast(c_ptr,1,MPI_FLOAT,0,MPI_COMM_WORLD);

}
