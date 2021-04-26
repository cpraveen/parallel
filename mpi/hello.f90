program mpitest
   use MPI
   integer :: rank, size, ierr, plen
   character*(MPI_MAX_PROCESSOR_NAME) :: pname

   call MPI_Init(ierr)
   call MPI_Comm_size(MPI_COMM_WORLD, size, ierr)
   call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
   call MPI_Get_processor_name(pname, plen, ierr)

   write(*,'("Hello World, from rank", I5, " of ", I5, " procs")', &
         ADVANCE="NO") rank, size
   write(*,'(", host ", A)') trim(pname)
   call MPI_Finalize(ierr)
end
