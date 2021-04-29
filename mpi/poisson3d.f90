!------------------------------------------------------------------------------
! 3d Poisson equation
!------------------------------------------------------------------------------
program main
   use MPI
   implicit none
   logical, dimension(1:3) :: pbc_check
   integer, dimension(1:3) :: spat_dim, proc_dim
   integer, dimension(1:3) :: loca_dim, mycoord
   integer, dimension(1:3) :: totmsgsize
   integer :: i, myid, numprocs, ierr, itermax, tag
   logical :: l_reorder
   integer :: req(1), status(MPI_STATUS_SIZE), GRID_COMM_WORLD
   integer :: myid_grid, nump_grid, tmp, t0, t1
   integer :: iStart, jStart, kStart, iEnd, jEnd, kEnd, MaxBufLen
   integer :: source, dest, dir, disp, iter
   double precision,parameter :: eps = 1.0d-10
   double precision :: maxdelta
   double precision,allocatable :: phi(:,:,:,:), fieldSend(:), fieldRecv(:)

   call MPI_Init(ierr)
   call MPI_Comm_rank(MPI_COMM_WORLD, myid, ierr)
   call MPI_Comm_size(MPI_COMM_WORLD, numprocs, ierr)
   if(myid.eq.0) then
      write(*,*) ' spat_dim , proc_dim, PBC ? '
      do i=1,3
         read(*,*) spat_dim(i), proc_dim(i), pbc_check(i)
      enddo
   endif

   call MPI_Bcast(spat_dim , 3, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
   call MPI_Bcast(proc_dim , 3, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
   call MPI_Bcast(pbc_check, 3, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)

   call MPI_Dims_create(numprocs, 3, proc_dim, ierr)
   if(myid.eq.0) write(*,'(a,3(i3,x))') 'Grid: ', (proc_dim(i),i=1,3)

   l_reorder = .true.
   call MPI_Cart_create(MPI_COMM_WORLD, 3, proc_dim, pbc_check, &
                        l_reorder, GRID_COMM_WORLD, ierr)

   if(GRID_COMM_WORLD .eq. MPI_COMM_NULL) goto 9
   call MPI_Comm_rank(GRID_COMM_WORLD, myid_grid, ierr)
   call MPI_Comm_size(GRID_COMM_WORLD, nump_grid, ierr)

   call MPI_Cart_coords(GRID_COMM_WORLD, myid_grid, 3, mycoord, ierr)

   do i=1,3
      loca_dim(i) = spat_dim(i)/proc_dim(i)
      if(mycoord(i) < mod(spat_dim(i),proc_dim(i))) then
         loca_dim(i) = loca_dim(i)+1
      endif
   enddo

   ! Solution variables
   iStart = 0 ; iEnd = loca_dim(3)+1
   jStart = 0 ; jEnd = loca_dim(2)+1
   kStart = 0 ; kEnd = loca_dim(1)+1
   allocate(phi(iStart:iEnd, jStart:jEnd, kStart:kEnd, 0:1))

   ! Buffers for send/recv
   ! j-k plane
   totmsgsize(3) = loca_dim(1)*loca_dim(2)
   MaxBufLen = max(MaxBufLen,totmsgsize(3))
   ! i-k plane
   totmsgsize(2) = loca_dim(1)*loca_dim(3)
   MaxBufLen = max(MaxBufLen,totmsgsize(2))
   ! i-j plane
   totmsgsize(1) = loca_dim(2)*loca_dim(3)
   MaxBufLen = max(MaxBufLen,totmsgsize(1))

   allocate(fieldSend(1:MaxBufLen))
   allocate(fieldRecv(1:MaxBufLen))

   ! Start iterations
   itermax = 100000
   t0 = 0 ; t1 = 1
   tag = 0
   do iter = 1, ITERMAX
      do disp = -1, 1, 2
         do dir = 1, 3
            call MPI_Cart_shift(GRID_COMM_WORLD, (dir-1), &
                                disp, source, dest, ierr)
            if(source /= MPI_PROC_NULL) then
               call MPI_Irecv(fieldRecv(1), totmsgsize(dir), &
                              MPI_DOUBLE_PRECISION, source, &
                              tag, GRID_COMM_WORLD, req(1), ierr)
            endif ! source exists

            if(dest /= MPI_PROC_NULL) then
               call CopySendBuf(phi(iStart, jStart, kStart, t0), &
                                iStart, iEnd, jStart, jEnd, kStart, kEnd, &
                                disp, dir, fieldSend, MaxBufLen)
               call MPI_Send(fieldSend(1), totmsgsize(dir), &
                             MPI_DOUBLE_PRECISION, dest, tag, &
                             GRID_COMM_WORLD, ierr)
            endif ! destination exists

            if(source /= MPI_PROC_NULL) then
               call MPI_Wait(req(1), status, ierr)
               call CopyRecvBuf(phi(iStart, jStart, kStart, t0), &
                                iStart, iEnd, jStart, jEnd, kStart, kEnd, &
                                disp, dir, fieldRecv, MaxBufLen)
            endif ! source exists
         enddo ! dir
      enddo   ! disp

      call Jacobi_sweep(loca_dim(1), loca_dim(2), loca_dim(3), & 
                        phi(iStart, jStart, kStart, 0), t0, t1, &
                        maxdelta)

      call MPI_Allreduce(MPI_IN_PLACE, maxdelta, 1, & 
                         MPI_DOUBLE_PRECISION, &
                         MPI_MAX, GRID_COMM_WORLD, ierr)
      if(maxdelta < eps) goto 9
      tmp=t0; t0=t1; t1=tmp
   enddo   ! iter

9  continue
   call MPI_Finalize(ierr)
end program main

!------------------------------------------------------------------------------
subroutine CopySendBuf(phi, iStart, iEnd, jStart, jEnd, kStart, kEnd, &
                       disp, dir, fieldSend, MaxBufLen)
   implicit none
   integer :: iStart, iEnd, jStart, jEnd, kStart, kEnd, disp, dir, MaxBufLen
   double precision :: phi(iStart:iEnd,jStart:jEnd,kStart:kEnd)
   double precision :: fieldSend(1:MaxBufLen)

end subroutine CopySendBuf

!------------------------------------------------------------------------------
subroutine CopyRecvBuf(phi, iStart, iEnd, jStart, jEnd, kStart, kEnd, &
                       disp, dir, fieldRecv, MaxBufLen)
   implicit none
   integer :: iStart, iEnd, jStart, jEnd, kStart, kEnd, disp, dir, MaxBufLen
   double precision :: phi(iStart:iEnd,jStart:jEnd,kStart:kEnd)
   double precision :: fieldRecv(1:MaxBufLen)
end subroutine CopyRecvBuf

!------------------------------------------------------------------------------
! Perform one iteration of Jacobi
!------------------------------------------------------------------------------
subroutine Jacobi_sweep(nx, ny, nz, phi, t0, t1, maxdelta)
   implicit none
   integer :: nx, ny, nz
   double precision :: phi(0:nx+1,0:ny+1,0:nz+1,0:2)
   integer :: t0, t1
   double precision :: maxdelta
   ! Local variables
   double precision,parameter :: const = 1.0d0/6.0d0
   double precision :: rhs, h ! TODO: h not defined
   integer :: i, j, k

   rhs = 1.0d0
   maxdelta = 0.0d0
   do k=1,nz
      do j=1,ny
         do i=1,nx
            phi(i,j,k,t1) = ( phi(i-1,j,k,t0) + phi(i+1,j,k,t0) &
                            + phi(i,j-1,k,t0) + phi(i,j+1,k,t0) &
                            + phi(i,j,k-1,t0) + phi(i,j,k+1,t0) &
                            + h**2 * rhs ) * const
         enddo
      enddo
   enddo
end subroutine Jacobi_sweep
