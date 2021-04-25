#include <stdio.h>
#include <mpi.h>

int main(int argc, char** argv) {
   int rank, size, ierror;

   ierror = MPI_Init(&argc, &argv);
   ierror = MPI_Comm_size(MPI_COMM_WORLD, &size);
   ierror = MPI_Comm_rank(MPI_COMM_WORLD, &rank);

   printf("Hello World, I am rank %d of %d procs\n", rank, size);
   ierror = MPI_Finalize();
   return 0;
}
