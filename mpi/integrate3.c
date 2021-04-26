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
{
   int size, rank, i, ierr;
   double a, b, res, mya, myb, psum;

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

   MPI_Reduce(&psum, &res, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
   if(rank == 0) printf("Result: %20.14e\n",res);

   MPI_Finalize();
}
