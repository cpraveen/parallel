c  Fortran example  
      program hello
      include 'mpif.h'
      integer rank, size, ierr
   
      call MPI_INIT(ierr)
      call MPI_COMM_SIZE(MPI_COMM_WORLD, size, ierr)
      call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)
      write(*,*) 'Hello World, I am rank ',rank,' of ',size,' procs'
      call MPI_FINALIZE(ierr)
      end
