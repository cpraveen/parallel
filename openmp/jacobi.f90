! Solve
!    -Laplace(u) = 1  in  (0,1)x(0,1)
!              u = 0  on  boundary
program main
   use omp_lib
   implicit none
   integer,parameter :: N = 500
   double precision, dimension(1:N,1:N,0:1) :: phi
   double precision, dimension(1:N,1:N) :: rhs
   double precision :: maxdelta, eps, h
   integer :: i, k, t0,t1, iter

!$OMP PARALLEL
   if(omp_get_thread_num() == 0) then
      print*,'Number of threads = ', omp_get_num_threads()
   endif
!$OMP END PARALLEL

   eps = 1.0d-10 ! convergence threshold
   h   = 1.0 / (N-1)
   t0 = 0; t1 = 1; iter = 0
   maxdelta = 2.0d0*eps

   phi = 0.0d0  ! Initial guess
   rhs = 1.0d0  ! rhs function

   do while(maxdelta.gt.eps)
      maxdelta = 0.d0
!$OMP PARALLEL DO COLLAPSE(2) REDUCTION(max:maxdelta)
      do k=2,N-1
         do i=2,N-1
            ! four flops, one store, four loads
            phi(i,k,t1) = ( phi(i+1,k,t0) + phi(i-1,k,t0) &
                          + phi(i,k+1,t0) + phi(i,k-1,t0) &
                          + h**2 * rhs(i,k)) * 0.25
            maxdelta = max(maxdelta, abs(phi(i,k,t1)-phi(i,k,t0)))
         enddo
      enddo
!$OMP END PARALLEL DO
      ! swap arrays
      i = t0 ; t0 = t1 ; t1 = i
      iter = iter + 1
      write(*,*) iter, maxdelta
   enddo

end program main
