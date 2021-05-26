#include <petsc.h>

int main(int argc, char **argv)
{
   PetscErrorCode ierr;
   PetscMPIInt    rank, size;
   ierr = PetscInitialize(&argc,&argv,NULL,"Hello World!"); CHKERRQ(ierr);
   ierr = MPI_Comm_rank(PETSC_COMM_WORLD,&rank); CHKERRQ(ierr);
   ierr = MPI_Comm_size(PETSC_COMM_WORLD,&size); CHKERRQ(ierr);
   ierr = PetscPrintf(PETSC_COMM_SELF, "Hello from rank %d of %d\n",rank,size); CHKERRQ(ierr);
   return PetscFinalize();
}
