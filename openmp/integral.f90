program integral
   use omp_lib
   double precision :: pi,w,sum,x
   integer :: i,n=1000000
   pi = 0.0d0
   w = 1.0d0/N
   sum = 0.0d0
   !$OMP PARALLEL PRIVATE(x) FIRSTPRIVATE(sum)
   if (omp_get_thread_num() == 0)then
      print*,'Number of threads = ', omp_get_num_threads()
   endif
   !$OMP DO
   do i=1,n
      x = w*(i-0.5d0)
      sum = sum + 4.0d0/(1.0d0+x*x)
   enddo
   !$OMP END DO
   !$OMP CRITICAL
   pi= pi + w*sum
   !$OMP END CRITICAL
   !$OMP END PARALLEL
   print*,'pi    = ', pi
   print*,'Error = ', abs(pi - 4.0d0*atan(1.0d0))
end program integral
