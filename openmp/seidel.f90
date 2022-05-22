! Solve
!    -Laplace(u) = 1  in  (0,1)x(0,1)
!              u = 0  on  boundary
! using Gauss-Seidel. See
!    Hager and Wallein, Intro to HPC for Scientists and Engineers, Section 6.3
program main
   use omp_lib
   implicit none
   integer,parameter :: N = 100
   integer,parameter :: imax = N, jmax = N, kmax = N
   double precision, allocatable :: phi(:,:,:), rhs(:,:,:)
   double precision,parameter :: osth = 1.0d0/6.0d0
   double precision,parameter :: eps = 1.0d-8
   double precision :: maxdelta, h, tmp
   integer :: threadId, numThreads, jStart, jEnd, i, j, k, l, iter

!$OMP PARALLEL
   numThreads = OMP_GET_NUM_THREADS()
   if(omp_get_thread_num() == 0) then
      print*,'Number of threads = ', numThreads
      if(mod(jmax, numThreads) /= 0)then
         print*, "jmax must be a multiple of numThreads"
      endif
   endif
!$OMP END PARALLEL

   if(mod(jmax, numThreads) /= 0)then
      print*,'jmax =', jmax
      stop
   endif

   allocate(phi(0:imax+1,0:jmax+1,0:kmax+1))
   allocate(rhs(0:imax+1,0:jmax+1,0:kmax+1))

   h   = 1.0 / (N-1)
   iter = 0
   maxdelta = 2.0d0*eps

   phi = 0.0d0  ! Initial guess
   rhs = 1.0d0  ! rhs function

   do while(maxdelta > eps)
      maxdelta = 0.d0

!$OMP PARALLEL PRIVATE(k,j,i,jStart,jEnd,threadID,tmp) reduction(max:maxdelta)
      threadID = OMP_GET_THREAD_NUM()
      jStart = jmax/numThreads*threadID
      jEnd = jStart + jmax/numThreads ! jmax is a multiple of numThreads
      do l=1,kmax+numThreads-1
         k = l - threadID
         if((k.ge.1).and.(k.le.kmax)) then
            do j=jStart,jEnd ! this is the actual parallel loop
               do i=1,imax
                  tmp = phi(i,j,k)
                  phi(i,j,k) = ( phi(i-1,j,k) + phi(i+1,j,k) &
                               + phi(i,j-1,k) + phi(i,j+1,k) &
                               + phi(i,j,k-1) + phi(i,j,k+1) &
                               + h**2 * rhs(i,j,k) ) * osth
               maxdelta = max(maxdelta, abs(tmp - phi(i,j,k)))
               enddo
            enddo
         endif
!$OMP BARRIER
      enddo
!$OMP END PARALLEL

      iter = iter + 1
      write(*,*) iter, maxdelta
   enddo

end program main
