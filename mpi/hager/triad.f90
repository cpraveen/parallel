!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! triad.f90 benchmark demo code
! G. Hager, 2010
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program triad
  implicit none

  double precision, dimension(:),allocatable :: A,B,C,D
! Intel-specific: 16-byte alignment of allocatables
!DEC$ ATTRIBUTES ALIGN: 16 :: A
!DEC$ ATTRIBUTES ALIGN: 16 :: B
!DEC$ ATTRIBUTES ALIGN: 16 :: C
!DEC$ ATTRIBUTES ALIGN: 16 :: D
  double precision :: MFLOPS,WT
  integer :: N,i
  integer(kind=8) :: R

  read *,N

  allocate(A(1:N),B(1:N),C(1:N),D(1:N))

  do i=1,N
     A(i) = 0.d0; B(i) = 1.d0
     C(i) = 2.d0; D(i) = 3.d0
  enddo

  R = 1

  ! warm up
  call do_triad(A,B,C,D,N,R,WT)

  do
     call do_triad(A,B,C,D,N,R,WT)
     ! exit if duration was above some limit
     if(WT.ge.0.2d0) exit
     ! else do it again with doubled repeat count
     R = R*2
  enddo

  MFLOPS = R*N*2.d0/(WT*1.d6) ! compute MFlop/sec rate
  print *, "Length: ",N,"   MFLOP/s: ",MFLOPS
  deallocate(A,B,C,D)
end program triad

subroutine do_triad(A,B,C,D,N,R,WT)
  implicit none
  integer, intent(in) :: N
  integer(kind=8), intent(in) :: R
  double precision, dimension(N), intent(out) :: A
  double precision, dimension(N), intent(in) :: B,C,D
  double precision, intent(out) :: WT
  double precision :: S,E
  integer :: N2
  ! assume 4MB outer level cache
  integer, parameter :: CACHE_LIMIT=131072
  integer :: i
  integer(kind=8) :: j

  N2 = N/2

  call get_walltime(S)
  
  if(N.le.CACHE_LIMIT) then
     do j=1,R
! Intel-specific: Assume aligned moves
!DEC$ vector aligned
!DEC$ vector temporal
        do i=1,N
           A(i) = B(i) + C(i) * D(i)
        enddo
        ! prevent loop interchange
        if(A(N2).lt.0) call dummy(A,B,C,D)
     enddo
  else
     do j=1,R
! Intel-specific: Assume aligned moves
!DEC$ vector aligned
!DEC$ vector nontemporal
        do i=1,N
           A(i) = B(i) + C(i) * D(i)
        enddo
        ! prevent loop interchange
        if(A(N2).lt.0) call dummy(A,B,C,D)
     enddo
  endif
  
  call get_walltime(E)
  
  WT = E-S

end subroutine do_triad

