program norm
   use omp_lib
   integer,parameter :: n=100
   integer :: i
   double precision :: a(n),  s

   do i=1,n
      a(i) = 1.0d0 * i
   enddo

   !$OMP PARALLEL DO REDUCTION(+:s)
   do i=1,n
      s = s + a(i)**2
   enddo
   !$OMP END PARALLEL DO
   print*,'norm    = ', dsqrt(s)
end program norm
