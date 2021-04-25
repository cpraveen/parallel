program main
   use MPI
   implicit none
   integer :: size, rank, i, ierror
   double precision :: a, b, res, mya, myb, psum, tmp
   double precision, parameter :: pi = 4.0d0 * atan(1.0d0)
   integer, dimension(MPI_STATUS_SIZE) :: status
   double precision, external :: integrate

   call MPI_Init(ierror)
   call MPI_Comm_size(MPI_COMM_WORLD, size, ierror)
   call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierror)

   ! integration limits
   a = 0.0d0 ; b = 2.0d0 * pi ; res = 0.0d0

   ! limits for "me"
   mya = a + rank*(b-a)/size
   myb = mya + (b-a)/size
   ! integrate f(x) over my own chunk - actual work
   psum = integrate(mya,myb)

   ! rank 0 collects partial results
   if(rank.eq.0) then
      print*,'Number of procs =', size
      res=psum
      do i=1,size-1
         call MPI_Recv(tmp, &
                       1,&
                       MPI_DOUBLE_PRECISION,&
                       i, &                   ! rank of source
                       0, &                   ! tag (unused here)
                       MPI_COMM_WORLD,&       ! communicator
                       status,&               ! status array (msg info)
                       ierror)
         res=res+tmp
      enddo
      write(*,*) 'Result: ',res
   ! ranks != 0 send their results to rank 0
   else
      call MPI_Send(psum, &                ! send buffer
                    1, &                   ! message length
                    MPI_DOUBLE_PRECISION,&
                    0, &                   ! rank of destination
                    0, &                   ! tag (unused here)
                    MPI_COMM_WORLD,ierror)
   endif

   call MPI_Finalize(ierror)
end program main

! Return integral of cos(x) from x=a to x=b
double precision function integrate(a, b)
   implicit none
   double precision :: a, b
   integrate = sin(b) - sin(a)
end function integrate
