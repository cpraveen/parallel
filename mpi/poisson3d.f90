!------------------------------------------------------------------------------
! 3d Poisson equation
!     -Laplace(u) = 1 in [0,1] x [0,1] x [0,1]
!              u  = 0 on boundary
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
   integer :: source, dest, dir, disp, iter, udim(2,3)
   double precision :: eps, maxdelta, h
   double precision,allocatable :: phi(:,:,:,:), fieldSend(:), fieldRecv(:)

   call MPI_Init(ierr)
   call MPI_Comm_rank(MPI_COMM_WORLD, myid, ierr)
   call MPI_Comm_size(MPI_COMM_WORLD, numprocs, ierr)
   ! Read parameters from file on rank=0
   if(myid.eq.0) then
      open(10,file='poisson3d.in',status='old')
      read(10,*) tmp
      read(10,*) proc_dim(1), proc_dim(2), proc_dim(3)
      read(10,*) itermax
      read(10,*) eps
      close(10)
      if(numprocs .ne. proc_dim(1)*proc_dim(2)*proc_dim(3))then
         print*,'Total procs cannot to factorized'
         print*,'Total procs = ', numprocs
         print*,'Proc grid   = ', proc_dim(:)
         call MPI_Abort(MPI_COMM_WORLD, tmp, ierr)
      endif
      spat_dim = [tmp, tmp, tmp]
      pbc_check = [.false., .false., .false.] ! Dirichlet bc
   endif

   ! Send parameters to all ranks
   call MPI_Bcast(spat_dim , 3, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
   call MPI_Bcast(proc_dim , 3, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
   call MPI_Bcast(pbc_check, 3, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
   call MPI_Bcast(itermax  , 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
   call MPI_Bcast(eps      , 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

   call MPI_Dims_create(numprocs, 3, proc_dim, ierr)

   ! We assume dx = dy = dz
   h = 1.0d0 / (spat_dim(1)-1)

   if(myid.eq.0)then
       write(*,'(a,3(i3,x))') 'Spatial Grid: ', (spat_dim(i),i=1,3)
       write(*,'(a,3(i3,x))') 'MPI     Grid: ', (proc_dim(i),i=1,3)
       write(*,'(a,e12.4)')   'Spatial h   : ', h
       write(*,'(a,i6)')      'itermax     : ', itermax
       write(*,'(a,e12.4)')   'eps         : ', eps
   endif

   l_reorder = .true.
   call MPI_Cart_create(MPI_COMM_WORLD, 3, proc_dim, pbc_check, &
                        l_reorder, GRID_COMM_WORLD, ierr)

   if(GRID_COMM_WORLD .eq. MPI_COMM_NULL)then
      if(myid.eq.0) print*,'Failed to create GRID_COMM_WORLD'
      call MPI_Abort(MPI_COMM_WORLD, tmp, ierr)
   endif

   call MPI_Comm_rank(GRID_COMM_WORLD, myid_grid, ierr)
   call MPI_Comm_size(GRID_COMM_WORLD, nump_grid, ierr)

   call MPI_Cart_coords(GRID_COMM_WORLD, myid_grid, 3, mycoord, ierr)

   ! loca_dim = grid size owned by current rank
   do i=1,3
      loca_dim(i) = spat_dim(i)/proc_dim(i)
      if(mycoord(i) < mod(spat_dim(i),proc_dim(i))) then
         loca_dim(i) = loca_dim(i)+1
      endif
   enddo

   ! Solution variables
   ! One layer of ghost points on all sides
   iStart = 0 ; iEnd = loca_dim(3)+1
   jStart = 0 ; jEnd = loca_dim(2)+1
   kStart = 0 ; kEnd = loca_dim(1)+1
   allocate(phi(iStart:iEnd, jStart:jEnd, kStart:kEnd, 0:1))

   ! Buffers for send/recv
   MaxBufLen = 0
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

   ! Find grid index range over which we must update solution
   disp = -1
   do dir = 1, 3
      call MPI_Cart_shift(GRID_COMM_WORLD, (dir-1), &
                          disp, source, dest, ierr)
      if(dest /= MPI_PROC_NULL) then ! we have a neighbour on left
         udim(1,dir) = 1
      else ! no neighbour on left, dirichlet bc
         udim(1,dir) = 2
      endif
      if(source /= MPI_PROC_NULL) then ! we have a neighbour on right
         udim(2,dir) = loca_dim(dir)
      else ! no neighbour on right, dirichlet bc
         udim(2,dir) = loca_dim(dir) - 1
      endif
   enddo

   ! Make initial guess
   ! Dirichlet bc is also zero
   phi = 0.0d0

   ! Start iterations
   maxdelta = 2.0d0 * eps
   t0 = 0 ; t1 = 1
   tag = 0
   iter = 0
   do while(iter < ITERMAX .and. maxdelta > eps)
      ! Each rank has at most two neighbours in every direction.
      ! Each rank sends data to one neighbour and receives from the other.
      !
      !     |-----|         |-----|         |-----|
      !  <--|  B  |<--   <--|  A  |<--   <--|  C  |<--
      !     |-----|         |-----|         |-----|
      !
      ! A sends to B and receives from C. At the same time B is waiting to 
      ! receive from A, and C is sending to A. Thus communication happens in
      ! correct order.
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

      call Jacobi_sweep(loca_dim(3), loca_dim(2), loca_dim(1), & 
                        phi(iStart, jStart, kStart, 0), t0, t1, &
                        udim, h, maxdelta)

      call MPI_Allreduce(MPI_IN_PLACE, maxdelta, 1, & 
                         MPI_DOUBLE_PRECISION, &
                         MPI_MAX, GRID_COMM_WORLD, ierr)
      iter = iter + 1
      if(myid.eq.0) print*,iter,maxdelta
      tmp=t0; t0=t1; t1=tmp ! swap solution values
   enddo  ! iter

   deallocate(phi, fieldSend, fieldRecv)
   call MPI_Finalize(ierr)
end program main

!------------------------------------------------------------------------------
! One layer of halo cells
!------------------------------------------------------------------------------
subroutine CopySendBuf(phi, iStart, iEnd, jStart, jEnd, kStart, kEnd, &
                       disp, dir, fieldSend, MaxBufLen)
   implicit none
   integer :: iStart, iEnd, jStart, jEnd, kStart, kEnd, disp, dir, MaxBufLen
   double precision :: phi(iStart:iEnd,jStart:jEnd,kStart:kEnd)
   double precision :: fieldSend(1:MaxBufLen)
   ! Local variables
   integer :: i1, i2, j1, j2, k1, k2, c
   integer :: i, j, k

   if(dir < 1 .or. dir > 3) stop 'CopySendBuf: dir is wrong'
   if(disp /= -1 .and. disp /= 1) stop 'CopySendBuf: disp is wrong'

   if(dir == 1)then ! i-j plane
      i1 = iStart+1; i2 = iEnd-1
      j1 = jStart+1; j2 = jEnd-1
      if(disp == -1)then
         k1 = 1; k2 = 1
      else
         k1 = kEnd-1; k2 = kEnd-1
      endif
   elseif(dir == 2)then ! i-k plane
      i1 = iStart+1; i2 = iEnd-1
      k1 = kStart+1; k2 = kEnd-1
      if(disp == -1)then
         j1 = 1; j2 = 1
      else
         j1 = jEnd-1; j2 = jEnd-1
      endif
   elseif(dir == 3)then ! j-k plane
      j1 = jStart+1; j2 = jEnd-1
      k1 = kStart+1; k2 = kEnd-1
      if(disp == -1)then
         i1 = 1; i2 = 1
      else
         i1 = iEnd-1; i2 = iEnd-1
      endif
   endif

   c = 1
   do k=k1,k2
      do j=j1,j2
         do i=i1,i2
            fieldSend(c) = phi(i,j,k)
            c = c + 1
         enddo
      enddo
   enddo

end subroutine CopySendBuf

!------------------------------------------------------------------------------
! One layer of halo cells
!------------------------------------------------------------------------------
subroutine CopyRecvBuf(phi, iStart, iEnd, jStart, jEnd, kStart, kEnd, &
                       disp, dir, fieldRecv, MaxBufLen)
   implicit none
   integer :: iStart, iEnd, jStart, jEnd, kStart, kEnd, disp, dir, MaxBufLen
   double precision :: phi(iStart:iEnd,jStart:jEnd,kStart:kEnd)
   double precision :: fieldRecv(1:MaxBufLen)
   ! Local variables
   integer :: i1, i2, j1, j2, k1, k2, c
   integer :: i, j, k

   if(dir < 1 .or. dir > 3) stop 'CopyRecvBuf: dir is wrong'
   if(disp /= -1 .and. disp /= 1) stop 'CopyRecvBuf: disp is wrong'

   if(dir == 1)then ! i-j plane
      i1 = iStart+1; i2 = iEnd-1
      j1 = jStart+1; j2 = jEnd-1
      if(disp ==  1)then
         k1 = 0; k2 = 0
      else
         k1 = kEnd; k2 = kEnd
      endif
   elseif(dir == 2)then ! i-k plane
      i1 = iStart+1; i2 = iEnd-1
      k1 = kStart+1; k2 = kEnd-1
      if(disp ==  1)then
         j1 = 0; j2 = 0
      else
         j1 = jEnd; j2 = jEnd
      endif
   elseif(dir == 3)then ! j-k plane
      j1 = jStart+1; j2 = jEnd-1
      k1 = kStart+1; k2 = kEnd-1
      if(disp ==  1)then
         i1 = 0; i2 = 0
      else
         i1 = iEnd; i2 = iEnd
      endif
   endif

   c = 1
   do k=k1,k2
      do j=j1,j2
         do i=i1,i2
            phi(i,j,k) = fieldRecv(c)
            c = c + 1
         enddo
      enddo
   enddo

end subroutine CopyRecvBuf

!------------------------------------------------------------------------------
! Perform one iteration of Jacobi
!------------------------------------------------------------------------------
subroutine Jacobi_sweep(nx, ny, nz, phi, t0, t1, udim, h, maxdelta)
   implicit none
   integer :: nx, ny, nz, udim(2,3)
   double precision :: phi(0:nx+1,0:ny+1,0:nz+1,0:1)
   integer :: t0, t1
   double precision :: h, maxdelta
   ! Local variables
   double precision,parameter :: const = 1.0d0/6.0d0
   double precision :: rhs
   integer :: i, j, k

   rhs = 1.0d0
   maxdelta = 0.0d0
   do k=udim(1,1),udim(2,1)
      do j=udim(1,2),udim(2,2)
         do i=udim(1,3),udim(2,3)
            phi(i,j,k,t1) = ( phi(i-1,j,k,t0) + phi(i+1,j,k,t0) &
                            + phi(i,j-1,k,t0) + phi(i,j+1,k,t0) &
                            + phi(i,j,k-1,t0) + phi(i,j,k+1,t0) &
                            + h**2 * rhs ) * const
            maxdelta = max(maxdelta, abs(phi(i,j,k,t1)-phi(i,j,k,t0)))
         enddo
      enddo
   enddo
end subroutine Jacobi_sweep
