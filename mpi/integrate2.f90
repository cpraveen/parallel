program main
   use MPI
   implicit none
   integer :: size, rank, i, ierror
   double precision :: a, b, res, mya, myb, psum
   double precision, parameter :: pi = 4.0d0 * atan(1.0d0)
   integer, allocatable, dimension(:,:) :: statuses
   integer, allocatable, dimension(:) :: requests
   double precision, allocatable, dimension(:) :: tmp
   double precision, external :: integrate

   call MPI_Init(ierror)
   call MPI_Comm_size(MPI_COMM_WORLD, size, ierror)
   call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierror)

   ! integration limits
   a = 0.0d0 ; b = 2.0d0 * pi ; res = 0.0d0

   if(rank.eq.0) then
      allocate(statuses(MPI_STATUS_SIZE, size-1))
      allocate(requests(size-1))
      allocate(tmp(size-1))
      ! pre-post nonblocking receives
      do i=1,size-1
         call MPI_Irecv(tmp(i), 1, MPI_DOUBLE_PRECISION, &
                        i, 0, MPI_COMM_WORLD,       &
                        requests(i), ierror)
      enddo
   endif
   ! limits for "me"
   mya = a + rank*(b-a)/size
   myb = mya + (b-a)/size

   ! integrate f(x) over my own chunk - actual work
   psum = integrate(mya,myb)

   ! rank 0 collects partial results
   if(rank.eq.0) then
      res = psum
      call MPI_Waitall(size-1, requests, statuses, ierror)
      do i=1,size-1
         res = res + tmp(i)
      enddo
      write (*,*) 'Result: ',res
   ! ranks != 0 send their results to rank 0
   else
      call MPI_Send(psum, 1, &
                    MPI_DOUBLE_PRECISION, 0, 0, &
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
