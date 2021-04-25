#include <stdio.h>
#include <mpi.h>

int main(int argc, char** argv) {
   int rank, size, ierr;

   ierr = MPI_Init(&argc, &argv);
   ierr = MPI_Comm_size(MPI_COMM_WORLD, &size);
   ierr = MPI_Comm_rank(MPI_COMM_WORLD, &rank);

   printf("Hello World, I am rank %d of %d procs\n", rank, size);
   ierr = MPI_Finalize();
   return 0;
}
