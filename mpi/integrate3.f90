program main
   use MPI
   implicit none
   integer :: size, rank, i, ierr
   double precision :: a, b, res, mya, myb, psum
   double precision, external :: integrate

   call MPI_Init(ierr)
   call MPI_Comm_size(MPI_COMM_WORLD, size, ierr)
   call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)

   ! integration limits
   a = 0.0d0 ; b = 2.0d0 ; res = 0.0d0

   ! limits for "me"
   mya = a + rank*(b-a)/size
   myb = mya + (b-a)/size

   ! integrate f(x) over my own chunk - actual work
   psum = integrate(mya,myb)

   call MPI_Reduce(psum, res, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, &
                   MPI_COMM_WORLD, ierr)

   if(rank.eq.0)then
      write(*,*)'Result =', res
   endif

   call MPI_Finalize(ierr)
end program main

! Return integral of cos(x) from x=a to x=b
function integrate(a, b)
   implicit none
   double precision :: a, b, integrate
   integrate = sin(b) - sin(a)
end function integrate
