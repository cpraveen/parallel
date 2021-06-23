/* Solves linear advection equation
             u_t + a u_x = 0
   with periodic bc and using MPI
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h>

double initial_condition(double x)
{
   return sin(8.0 * M_PI * x);
}

void get_ghost_values(double size, double rank, int n1, double* u)
{
   int left, right;
   MPI_Request request[4];
   MPI_Status status[4];

   if(size == 1) // serial run: apply periodicity
   {
      u[0]    = u[n1];
      u[n1+1] = u[1];
      return;
   }

   // Determine left/right ranks by periodicity
   if(rank == 0)
   {
      left  = size-1;
      right = rank + 1;
   }
   else if(rank == size-1)
   {
      left  = rank-1;
      right = 0;
   }
   else
   {
      left  = rank - 1;
      right = rank + 1;
   }

   MPI_Isend(&u[1],1,MPI_DOUBLE,left,0,MPI_COMM_WORLD,&request[0]);
   MPI_Isend(&u[n1],1,MPI_DOUBLE,right,0,MPI_COMM_WORLD,&request[1]);
   MPI_Irecv(&u[0],1,MPI_DOUBLE,left,0,MPI_COMM_WORLD,&request[2]);
   MPI_Irecv(&u[n1+1],1,MPI_DOUBLE,right,0,MPI_COMM_WORLD,&request[3]);
   MPI_Waitall(4, request, status);
}

void save_solution(int rank, int n1, double* x, double* u)
{
   static int count = 0;
   int i;
   FILE *fp;
   char fname[80];
   if(count == 0)
      sprintf(fname, "ini_%03d", rank);
   else
      sprintf(fname, "sol_%03d", rank);
   strcat(fname, ".txt");

   fp = fopen(fname, "w");
   for(i=1; i<=n1; ++i)
   {
      fprintf(fp, "%e %e\n", x[i], u[i]);
   }
   fclose(fp);
   ++count;
}

// Grid of n points is constructed in the domain [xmin,xmax] such that
// x[0]   = xmin
// x[n-1] = xmax - dx
int main(int argc, char** argv)
{
   double a = 1.0; // advection speed
   int n = 100;    // number of grid points
   int n1, size, rank, i, ierr, iter;
   double xmin, xmax, xmin1, xmax1, dx, cfl=0.95, t, Tf=2.0, dt, nu;
   double *x, *u, *unew;

   MPI_Init(&argc, &argv);
   MPI_Comm_size(MPI_COMM_WORLD, &size);
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);

   xmin = 0.0; xmax = 1.0;
   dx = (xmax - xmin)/n;
   dt = cfl * fabs(a) * dx;

   // find local grid size and local domain
   n1 = n / size;
   xmin1 = xmin + dx * rank * n1;
   if(rank == size-1) n1 = n - (size-1)*n1;
   xmax1 = xmin1 + dx * n1;

   x    = (double*)malloc((n1+2) * sizeof(double));
   u    = (double*)malloc((n1+2) * sizeof(double));
   unew = (double*)malloc((n1+2) * sizeof(double));

   // set initial condition
   for(i=1; i<=n1; ++i)
   {
      x[i] = xmin1 + (i-1) * dx;
      u[i] = initial_condition(x[i]);
   }

   // Get ghost values
   get_ghost_values(size, rank, n1, u);

   // save solution
   save_solution(rank, n1, x, u);

   t = 0.0; iter = 0;
   while(t < Tf)
   {
      if(t + dt > Tf) dt = Tf - t;
      nu = a*dt/dx;

      // update solution using Lax-Wendroff
      for(i=1; i<=n1; ++i)
      {
         unew[i] = u[i] - 0.5*nu*(u[i+1] - u[i-1]) 
                        + 0.5*nu*nu*(u[i-1] - 2.0*u[i] + u[i+1]);
      }

      // swap solution
      for(i=1; i<=n1; ++i) u[i] = unew[i];

      // Get ghost values
      get_ghost_values(size, rank, n1, u);

      t += dt; ++iter;
      if(rank == 0) printf("iter, t = %5d, %12.6e\n", iter, t);
   }

   // save solution
   save_solution(rank, n1, x, u);
   if(rank == 0) printf("Run: python linadv1d.py\n");

   free(u); free(unew);
   MPI_Finalize();
   return 0;
}
