program mpitest
   use MPI
   integer :: rank, size, ierr

   call MPI_Init(ierr)
   call MPI_Comm_size(MPI_COMM_WORLD, size, ierr)
   call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)

   write(*,*) 'Hello World, I am rank ',rank,' of ',size,' procs'
   call MPI_Finalize(ierr)
end
