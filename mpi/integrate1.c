#include <stdio.h>
#include <math.h>
#include <mpi.h>

// Return integral of cos(x) from x=a to x=b
double integrate(double a, double b)
{
   return sin(b) - sin(a);
}

// Compute integral(x=a to x=b) f(x)
int main(int argc, char** argv)
   int size, rank, i, ierr;
   double a, b, res, mya, myb, psum, tmp;
   MPI_Status status;

   MPI_Init(&argc, &argv);
   MPI_Comm_size(MPI_COMM_WORLD, &size);
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);

   // integration limits
   a = 0.0; b = 2.0; res = 0.0;

   // limits for "me"
   mya = a + rank*(b-a)/size;
   myb = mya + (b-a)/size;
   // integrate f(x) over my own chunk - actual work
   psum = integrate(mya, myb);

   // rank 0 collects partial results
   if(rank == 0)
   {
      printf("Number of procs = %d\n", size);
      res = psum;
      for(i=1; i<size; ++i)
      {
         MPI_Recv(&tmp,
                  1,
                  MPI_DOUBLE_PRECISION,
                  i,                    // rank of source
                  0,                    // tag (unused here)
                  MPI_COMM_WORLD,       // communicator
                  &status);             // status array (msg info)
         res += tmp;
      }
      printf("Result: %e\n",res);
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

   MPI_Finalize();
}
