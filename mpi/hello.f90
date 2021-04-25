program mpitest
   use MPI
   integer :: rank, size, ierror

   call MPI_Init(ierror)
   call MPI_Comm_size(MPI_COMM_WORLD, size, ierror)
   call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierror)

   write(*,*) 'Hello World, I am rank ',rank,' of ',size,' procs'
   call MPI_Finalize(ierror)
end
