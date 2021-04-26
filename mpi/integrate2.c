#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>

// Return integral of cos(x) from x=a to x=b
double integrate(double a, double b)
{
   return sin(b) - sin(a);
}

// Compute integral(x=a to x=b) f(x)
int main(int argc, char** argv)
{
   int size, rank, i, ierr;
   double a, b, res, mya, myb, psum, *tmp;
   MPI_Request *request;
   MPI_Status  *status;

   MPI_Init(&argc, &argv);
   MPI_Comm_size(MPI_COMM_WORLD, &size);
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);

   // integration limits
   a = 0.0; b = 2.0; res = 0.0;

   // Start non-blocking recives on rank=0
   if(rank == 0)
   {
      status  = (MPI_Status*)malloc(size*sizeof(MPI_Status));
      request = (MPI_Request*)malloc(size*sizeof(MPI_Request));
      tmp     = (double*)malloc(size*sizeof(double));
      for(i=1; i<size; ++i)
         MPI_Irecv(&tmp[i], 1, MPI_DOUBLE, i, 0, MPI_COMM_WORLD, &request[i]);
   }

   // limits for "me"
   mya = a + rank*(b-a)/size;
   myb = mya + (b-a)/size;
   // integrate f(x) over my own chunk - actual work
   psum = integrate(mya, myb);

   // rank 0 collects partial results
   if(rank == 0)
   {
      res = psum;
      MPI_Waitall(size-1, &request[1], &status[1]);
      for(i=1; i<size; ++i)
         res += tmp[i];
      printf("Result: %20.14e\n",res);
   }
   // ranks != 0 send their results to rank 0
   else
   {
      MPI_Send(&psum,                // send buffer
               1,                    // message length
               MPI_DOUBLE_PRECISION,
               0,                    // rank of destination
               0,                    // tag (unused here)
               MPI_COMM_WORLD);
   }

   if(rank == 0)
   {
      free(request); free(status); free(tmp);
   }
   MPI_Finalize();
}
