C
C     setenv OMP_SCHEDULE static ; setenv OMP_NUM_THREADS 3 ; env PINOMP_MASK=6 mpirun -np 2 
C     -npernode 1 /apps/rrze/bin/pin_omp -c "0,2,3" ./omp1 < input   
C
      PROGRAM JacobiMPI
      IMPLICIT NONE

      include 'mpif.h'
      integer myid, numprocs, ierr, myid_grid, nump_grid

      integer(KIND=4), allocatable, dimension(:) :: itest
      real(KIND=8), allocatable, dimension(:,:,:) :: fieldA, fieldB
      logical pbc_check(1:3), l_reorder

      integer spat_dim(1:3), proc_dim(1:3)
     $     , loca_dim(1:3), my_coords(1:3)

      integer n_dim , GRID_COMM_WORLD

      integer iStart, jStart, kStart, iEnd, jEnd, kEnd  

      integer tl_required, tl_provided

      integer i 
      real(KIND=8) tmp
C
C     Init MPI
C
C      call MPI_INIT( ierr )

C
C     Init HYBRID
C
      call MPI_INIT_THREAD( MPI_THREAD_MULTIPLE , tl_provided , ierr)

C
C     Get MPI Rank
C
      call MPI_COMM_RANK( MPI_COMM_WORLD, myid, ierr )      
C
C     How many processes are out there?
C
      call MPI_COMM_SIZE( MPI_COMM_WORLD, numprocs, ierr )  

      if(myid.eq.0) then
         write(*,'(a,i8,a,4(i8,x))') "# Threadlevel provided: ",tl_provided,
     $        " aus ",MPI_THREAD_SINGLE,MPI_THREAD_FUNNELED,
     $        MPI_THREAD_SERIALIZED, MPI_THREAD_MULTIPLE
         write(*,'(a)') "# Hybrid Version1: Master updates boundary cells and does halo exchange " 
C 2D does not really work
         n_dim=3
      endif

      call MPI_BCAST(n_dim , 1,MPI_INTEGER  ,0,MPI_COMM_WORLD,ierr)

      do i=1,n_dim
         spat_dim(i) = -1
         loca_dim(i) = -1
         proc_dim(i) = 0
      enddo

      if(myid.eq.0) then

         write(*,*) " spat_dim , proc_dim, PBC ? "
         do i=1,n_dim
            read(*,*) spat_dim(i) , proc_dim(i), pbc_check(i)
            write(*,*) i,"-Dim Input ",spat_dim(i) , proc_dim(i), pbc_check(i)
         enddo

c         pbc_check(1) = .false.
c         pbc_check(2) = .false.
c         pbc_check(3) = .false.

      endif

      call MPI_BCAST(spat_dim , 3,MPI_INTEGER  ,0,MPI_COMM_WORLD,ierr)

      call MPI_BCAST(proc_dim , 3,MPI_INTEGER  ,0,MPI_COMM_WORLD,ierr)

      call MPI_BCAST(pbc_check, 3,MPI_LOGICAL  ,0,MPI_COMM_WORLD,ierr)

      call MPI_DIMS_CREATE( numprocs , n_dim ,  proc_dim, ierr)

      if(myid.eq.0) write(*,'(a,i3,a,3(i3,x))') "# ",myid," Prozessorgrid ",
     $     (proc_dim(i),i=1,n_dim)


      l_reorder = .true.
      call MPI_CART_CREATE(MPI_COMM_WORLD, n_dim, proc_dim , pbc_check, 
     $     l_reorder, GRID_COMM_WORLD, ierr )


      if(GRID_COMM_WORLD .eq. MPI_COMM_NULL) goto 999

C      write(*,*),myid,GRID_COMM_WORLD, MPI_COMM_NULL

C
C     Get MPI Rank
C
      call MPI_COMM_RANK( GRID_COMM_WORLD, myid_grid, ierr )     
cc      write(42,*) myid, myid_grid,ierr
C
C     How many processes are still out there?
C
      call MPI_COMM_SIZE(GRID_COMM_WORLD , nump_grid, ierr )  
C
C     
      call MPI_Cart_coords(GRID_COMM_WORLD, myid_grid,3,
     $             my_coords,ierr)
      do i=1,n_dim
         loca_dim(i) = spat_dim(i) / proc_dim(i)
         if(my_coords(i) < mod(spat_dim(i),proc_dim(i))) then
            loca_dim(i) = loca_dim(i)+1
         endif
      enddo

      iStart = 0
      jStart = 0
      kStart = 0

      iEnd   = loca_dim(3)+1
      jEnd   = loca_dim(2)+1
      kEnd   = loca_dim(1)+1

C
C     Allocate fieldA including 1 boundary layer
C
      allocate(fieldA(iStart:iEnd ,jStart:jEnd,kStart:kEnd ))

C
C     Allocate fieldB including 1 boundary layer
C
      allocate(fieldB(iStart:iEnd ,jStart:jEnd,kStart:kEnd ))

      if(myid.eq.0) then
         tmp = 1.d0
         do i=1,n_dim
            tmp= tmp *dble(loca_dim(i)+2)
         enddo
         
         tmp= 2.d0*8.d0 * tmp /(1000.d0*1000.d0)

         write(*,'(a,i4,x,a,g12.6,x,a)') "#",myid,
     $       " Allocated 2 arrays with a total of " , tmp, " Mbytes "
      endif

      call ParallelSolver(n_dim, spat_dim, loca_dim, proc_dim, 
     $     fieldA , fieldB, iStart, iEnd, jStart, jEnd, kStart, kEnd, 
     $     GRID_COMM_WORLD, myid_grid, nump_grid)

 999  continue

      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
C
C     End MPI
C
      call MPI_FINALIZE(ierr)

      end


      subroutine ParallelSolver(n_dim, spat_dim, loca_dim, proc_dim, 
     $     fieldA , fieldB, iStart, iEnd, jStart, jEnd, kStart, kEnd, 
     $     GRID_COMM_WORLD, myid_grid, nump_grid)

      IMPLICIT NONE

      
      include 'omp_lib.h'
      include 'mpif.h'

      integer ierr, myid_grid, nump_grid , GRID_COMM_WORLD

      integer n_dim, spat_dim(1:3), proc_dim(1:3), loca_dim(1:3) , my_offset(1:3)

      integer iStart, jStart, kStart, iEnd, jEnd, kEnd  
      real(KIND=8) fieldA(iStart:iEnd ,jStart:jEnd,kStart:kEnd) , 
     $     fieldB(iStart:iEnd ,jStart:jEnd,kStart:kEnd)

      real(KIND=8), allocatable, dimension(:) :: fieldSend, fieldRecv
      integer MaxBufLen
     

      integer mycoord(1:3) ,i , j, k, counter , iter
      integer source, dest

      integer sendmsg(2,3,3,0:1),recvmsg(2,3,3,0:1)

      integer totmsgsize(1:3), dir, dir_mpi, tag , req(6)

      integer status(MPI_STATUS_SIZE)
      real*8 starttime, runtime,maxtime,mintime
      real*8 JacobiTime, JacobiTimeMi, JacobiTimeMa

      real*8 tmp, residual


      integer  ompid, numthreads
      real*8 JacobiTimeT
      integer iStartT, jStartT, kStartT, iEndT, jEndT, kEndT, workpT


      call InitBufferDims(1,loca_dim(1),1,loca_dim(2),1,loca_dim(3)
     $     ,totmsgsize,sendmsg,recvmsg,MaxBufLen)

C
C     Allocate fieldSend: SendBuffer holding NO boundary layer
C
      allocate(fieldSend(1:MaxBufLen))
C
C     Allocate fieldRecv: ReceiveBuffer holding NO boundary layer
C
      allocate(fieldRecv(1:MaxBufLen))
      
      call MPI_CART_COORDS(GRID_COMM_WORLD, myid_grid, 3 , mycoord ,ierr)

c      write(42+myid_grid,'(a,i4,a,3(i2,x),4x,3(i4,x))') "# ",myid_grid, 
c     $     " Coordinates ",
c     $     (mycoord(i),i=1,n_dim),kEnd,jEnd,iEnd
C
C     Offsets für globale Indizes
C
      do i=1,n_dim
         my_offset(i) = (spat_dim(i) / proc_dim(i)) * mycoord(i)
      enddo
c      write(42+myid_grid,'(a,i4,a,3(i6,x))') "# ",myid_grid, 
c     $     " Offsets ",
c     $     (my_offset(i),i=1,n_dim)
      
C
C     Initialize
C
      call InitArrays(loca_dim(1),loca_dim(2),loca_dim(3), 
     $     fieldA, fieldB, proc_dim , mycoord, my_offset)
C
C     
C
      do k = kStart+1, kEnd-1
         do j = jStart+1, jEnd-1
            do i = iStart+1, iEnd-1
               fieldA(i,j,k) = 0.d0
            enddo
         enddo
      enddo
C
C Residual    
C
      tmp=0.d0

      do k = kStart+1, kEnd-1
         do j = jStart+1, jEnd-1
            do i = iStart+1, iEnd-1
               tmp=tmp+(fieldA(i,j,k)-fieldB(i,j,k))*
     $              (fieldA(i,j,k)-fieldB(i,j,k))
            enddo
         enddo
      enddo

      call MPI_REDUCE(tmp , residual , 1 , MPI_DOUBLE_PRECISION, 
     $     MPI_SUM, 0, GRID_COMM_WORLD, ierr)

      if(myid_grid.eq.0) write(*,'(a,g20.12)') "# StartResidual ",residual 

C
C     Solve Equation (max 10 iterations)
C

      do iter = 1,10

         call MPI_BARRIER(GRID_COMM_WORLD , ierr)
         starttime=MPI_WTIME()

         ompid = 0
         numthreads = 1
         JacobiTime=0.d0

!$OMP PARALLEL private(ompid,JacobiTimeT,iStartT, jStartT, kStartT, iEndT, jEndT, kEndT,workpT)         

         JacobiTimeT=0.d0

         
!$       ompid       = omp_get_thread_num()
!$       numthreads  = omp_get_num_threads()
         if(ompid.ne.0) then
C
C     All Threads - except Master compute
C     
            kstartT=2
            kendT  =loca_dim(1)-1
C
C     Multiple Workers
C
            workpT = ((loca_dim(1)-1) - 2 +1) / (numthreads-1)

            kStartT=2 + (ompid-1) * workpT
            kEndT  =2 +  ompid    * workpT - 1

            if(ompid.eq.(numthreads-1)) kEndT = loca_dim(1)-1
C            write(*,*) " Start PARALLEL ",numthreads,ompid,kStartT, kEndT           


            jstartT=2
            jendT  =loca_dim(2)-1

            istartT=2
            iendT  =loca_dim(3)-1


            if(mod(iter,2).eq.1) then
               call Jacobi_solve1(loca_dim(1),loca_dim(2),loca_dim(3), fieldA, fieldB,
     $              kstartT,kendT,jstartT,jendT,istartT,iendT)
            else
               call Jacobi_solve1(loca_dim(1),loca_dim(2),loca_dim(3), fieldB, fieldA,
     $              kstartT,kendT,jstartT,jendT,istartT,iendT)
            endif

            JacobiTimeT=MPI_WTIME()-starttime

         else

            do dir = 0 , 1

               dir_mpi = 2*dir -1

               do counter = 1 , 3
               
                  call MPI_CART_SHIFT(GRID_COMM_WORLD , (counter-1) , dir_mpi , source, dest, ierr)
               
                  tag= dir_mpi * counter + 2*counter
                  if(source.ne.MPI_PROC_NULL) 
     $                 call MPI_IRECV(fieldRecv(1), totmsgsize(counter) , 
     $                 MPI_DOUBLE_PRECISION, source, 
     $                 tag ,GRID_COMM_WORLD, req(1), ierr)
                  
                  if(dest.ne.MPI_PROC_NULL) then
                     
                     if(mod(iter,2).eq.1) then
                        call CopySendBuf(fieldB,iStart,iEnd,jStart,jEnd,
     $                       kStart,kEnd, fieldSend,MaxBufLen, 
     $                       sendmsg(1,1,counter,dir),fieldA)
                     else
                        call CopySendBuf(fieldA,iStart,iEnd,jStart,jEnd,
     $                       kStart,kEnd, fieldSend,MaxBufLen, 
     $                       sendmsg(1,1,counter,dir),fieldB)
                     endif
                  
                     
                     call MPI_SEND(fieldSend(1), totmsgsize(counter) , MPI_DOUBLE_PRECISION, 
     $                    dest, tag, GRID_COMM_WORLD, ierr)
                     
                  endif

                  if(source.ne.MPI_PROC_NULL) then
                     call MPI_WAIT(req, status, ierr)
                     
                     if(mod(iter,2).eq.1) then
                        call CopyRecvBuf(fieldB,iStart,iEnd,jStart,jEnd,
     $                       kStart,kEnd, fieldRecv,MaxBufLen, 
     $                       recvmsg(1,1,counter,dir))
                     else
                        call CopyRecvBuf(fieldA,iStart,iEnd,jStart,jEnd,
     $                       kStart,kEnd, fieldRecv,MaxBufLen, 
     $                       recvmsg(1,1,counter,dir))
                     endif
                     
                  endif

               enddo

            enddo

         endif

!$OMP CRITICAL
         
         JacobiTime = max( JacobiTime , JacobiTimeT )

!$OMP END CRITICAL     

!$OMP END PARALLEL

         call MPI_BARRIER(GRID_COMM_WORLD , ierr)
         runtime=MPI_WTIME()-starttime
         call MPI_REDUCE(runtime , maxtime , 1 , MPI_DOUBLE_PRECISION, 
     $        MPI_MAX, 0, GRID_COMM_WORLD, ierr)
         call MPI_REDUCE(runtime , mintime , 1 , MPI_DOUBLE_PRECISION, 
     $        MPI_MIN, 0, GRID_COMM_WORLD, ierr)
         
         call MPI_REDUCE(JacobiTime , JacobiTimeMi , 1 , MPI_DOUBLE_PRECISION, 
     $        MPI_MIN, 0, GRID_COMM_WORLD, ierr)
         call MPI_REDUCE(JacobiTime , JacobiTimeMa , 1 , MPI_DOUBLE_PRECISION, 
     $        MPI_MAX, 0, GRID_COMM_WORLD, ierr)
         
C
C Residual    
C
         tmp=0.d0

         do k = kStart+1, kEnd-1
            do j = jStart+1, jEnd-1
               do i = iStart+1, iEnd-1
                  tmp=tmp+(fieldA(i,j,k)-fieldB(i,j,k))*
     $                 (fieldA(i,j,k)-fieldB(i,j,k))
               enddo
            enddo
         enddo
         
         call MPI_REDUCE(tmp , residual , 1 , MPI_DOUBLE_PRECISION, 
     $        MPI_SUM, 0, GRID_COMM_WORLD, ierr)
         
         tmp = 1.d0
         do i=1,n_dim
            tmp= tmp *dble(spat_dim(i))
         enddo
         
         tmp= tmp /(1000.d0*1000.d0*maxtime)
         
         if(myid_grid.eq.0) write(*,'(i3,x,a,4(g20.12,x),2(a,g20.12))') i,
     $        " Maxtime , Mintime + JacobiMi , JacobiMa ",maxtime,mintime,
     $        JacobiTimeMi,JacobiTimeMa,"  Residual ",residual, "  MLUPs ",tmp
         
      enddo


      
      return
      end


      subroutine InitBufferDims(kStartB,kendB,jStartB,jEndB,
     $     iStartB,iEndB,totmsgsize
     $     ,sendmsg,recvmsg,MaxBufLen)

      implicit none
      integer iStartB,iEndB,jStartB,jEndB,kStartB,kEndB
      integer totmsgsize(1:3)
      integer sendmsg(2,3,3,0:1),recvmsg(2,3,3,0:1)
      integer MaxBufLen

      MaxBufLen = -1
C
C     k-j plane
C
      totmsgsize(3) = (kEndB - kStartB+1)*  (jEndB - jStartB+1)
      MaxBufLen=max(MaxBufLen,totmsgsize(3))

      totmsgsize(2) = (kEndB - kStartB+1)*  (iEndB - iStartB+1)
      MaxBufLen=max(MaxBufLen,totmsgsize(2))

      totmsgsize(1) = (jEndB - jStartB+1)*  (iEndB - iStartB+1)
      MaxBufLen=max(MaxBufLen,totmsgsize(1))

      sendmsg(1,3,3,1) = iEndB
      sendmsg(2,3,3,1) = iEndB

      sendmsg(1,2,3,1) = 1
      sendmsg(2,2,3,1) = jEndB

      sendmsg(1,1,3,1) = 1
      sendmsg(2,1,3,1) = kEndB

      sendmsg(1,3,2,1) = 1
      sendmsg(2,3,2,1) = iEndB

      sendmsg(1,2,2,1) = jEndB
      sendmsg(2,2,2,1) = jEndB

      sendmsg(1,1,2,1) = 1
      sendmsg(2,1,2,1) = kEndB

      sendmsg(1,3,1,1) = 1
      sendmsg(2,3,1,1) = iEndB

      sendmsg(1,2,1,1) = 1
      sendmsg(2,2,1,1) = jEndB

      sendmsg(1,1,1,1) = kEndB
      sendmsg(2,1,1,1) = kEndB


      sendmsg(1,3,3,0) = 1
      sendmsg(2,3,3,0) = 1

      sendmsg(1,2,3,0) = 1
      sendmsg(2,2,3,0) = jEndB

      sendmsg(1,1,3,0) = 1
      sendmsg(2,1,3,0) = kEndB

      sendmsg(1,3,2,0) = 1
      sendmsg(2,3,2,0) = iEndB

      sendmsg(1,2,2,0) = 1
      sendmsg(2,2,2,0) = 1

      sendmsg(1,1,2,0) = 1
      sendmsg(2,1,2,0) = kEndB

      sendmsg(1,3,1,0) = 1
      sendmsg(2,3,1,0) = iEndB

      sendmsg(1,2,1,0) = 1
      sendmsg(2,2,1,0) = jEndB

      sendmsg(1,1,1,0) = 1
      sendmsg(2,1,1,0) = 1



      recvmsg(1,3,3,0) = iEndB+1
      recvmsg(2,3,3,0) = iEndB+1

      recvmsg(1,2,3,0) = 1
      recvmsg(2,2,3,0) = jEndB

      recvmsg(1,1,3,0) = 1
      recvmsg(2,1,3,0) = kEndB

      recvmsg(1,3,2,0) = 1
      recvmsg(2,3,2,0) = iEndB

      recvmsg(1,2,2,0) = jEndB+1
      recvmsg(2,2,2,0) = jEndB+1

      recvmsg(1,1,2,0) = 1
      recvmsg(2,1,2,0) = kEndB

      recvmsg(1,3,1,0) = 1
      recvmsg(2,3,1,0) = iEndB

      recvmsg(1,2,1,0) = 1
      recvmsg(2,2,1,0) = jEndB

      recvmsg(1,1,1,0) = kEndB+1
      recvmsg(2,1,1,0) = kEndB+1


      recvmsg(1,3,3,1) = 0
      recvmsg(2,3,3,1) = 0

      recvmsg(1,2,3,1) = 1
      recvmsg(2,2,3,1) = jEndB

      recvmsg(1,1,3,1) = 1
      recvmsg(2,1,3,1) = kEndB

      recvmsg(1,3,2,1) = 1
      recvmsg(2,3,2,1) = iEndB

      recvmsg(1,2,2,1) = 0
      recvmsg(2,2,2,1) = 0

      recvmsg(1,1,2,1) = 1
      recvmsg(2,1,2,1) = kEndB

      recvmsg(1,3,1,1) = 1
      recvmsg(2,3,1,1) = iEndB

      recvmsg(1,2,1,1) = 1
      recvmsg(2,2,1,1) = jEndB

      recvmsg(1,1,1,1) = 0
      recvmsg(2,1,1,1) = 0


      return
      end


      subroutine CopySendBuf(field,iStart,iEnd,jStart,jEnd,
     $     kStart,kEnd, fieldSend,MaxBufLen, 
     $     loopbonds,oldvals)

      implicit none

      integer iStart,iEnd,jStart,jEnd,kStart,kEnd
      real(KIND=8) field(iStart:iEnd ,jStart:jEnd,kStart:kEnd) 
      real(KIND=8) oldvals(iStart:iEnd ,jStart:jEnd,kStart:kEnd) 

      integer MaxBufLen
      real(KIND=8) fieldSend(1:MaxBufLen) 

      integer loopbonds(2,3)

      integer i,j,k
      integer offset

      do k=loopbonds(1,1),loopbonds(2,1)

         offset = 1-loopbonds(1,3) + (k-loopbonds(1,1))*
     $        (loopbonds(2,2)-loopbonds(1,2)+1)*
     $        (loopbonds(2,3)-loopbonds(1,3)+1)
            
         do j=loopbonds(1,2),loopbonds(2,2)
               
            do i=loopbonds(1,3),loopbonds(2,3)
               
               field(i,j,k)       = 0.1666666666d0* ( oldvals(i-1,j,k) + oldvals(i+1,j,k) 
     $              + oldvals(i,j-1,k) + oldvals(i,j+1,k) 
     $              + oldvals(i,j,k-1) + oldvals(i,j,k+1) )
               
               fieldSend(i+offset) = field(i,j,k)
               
            enddo

            offset=offset+loopbonds(2,3)-loopbonds(1,3)+1

         enddo


      enddo

      return
      end

      subroutine CopyRecvBuf(field,iStart,iEnd,jStart,jEnd,
     $     kStart,kEnd, fieldRecv,MaxBufLen, 
     $     loopbonds)

      implicit none

      integer iStart,iEnd,jStart,jEnd,kStart,kEnd
      real(KIND=8) field(iStart:iEnd ,jStart:jEnd,kStart:kEnd) 

      integer MaxBufLen
      real(KIND=8) fieldRecv(1:MaxBufLen) 

      integer loopbonds(2,3)

      integer i,j,k
      integer offset

      do k=loopbonds(1,1),loopbonds(2,1)

         offset = 1-loopbonds(1,3) + (k-loopbonds(1,1))*
     $        (loopbonds(2,2)-loopbonds(1,2)+1)*
     $        (loopbonds(2,3)-loopbonds(1,3)+1)
         
         do j=loopbonds(1,2),loopbonds(2,2)
            
            do i=loopbonds(1,3),loopbonds(2,3)
               
               field(i,j,k) = fieldRecv(i+offset)
               
            enddo
            
            offset=offset+loopbonds(2,3)-loopbonds(1,3)+1
            
         enddo

      enddo

      return
      end

      subroutine Jacobi_solve1(Nk,Nj,Ni, x, y,kstart,kend,jstart,jend,istart,iend)
      implicit none
      include 'omp_lib.h'

      integer kstart,kend,jstart,jend,istart,iend
      integer Nk,Nj,Ni
      real*8 x(0:Ni+1,0:Nj+1,0:Nk+1)
      real*8 y(0:Ni+1,0:Nj+1,0:Nk+1)

      integer i,j,k, jb, jblock

      jblock = 50

      do jb = jstart,jend,jblock
      
         do k = kstart,kend

            do j = jb , min(jb+jblock-1,jend)

               do i = istart,iend

                  y(i,j,k) = 0.1666666666d0* ( x(i-1,j,k) + x(i+1,j,k) 
     $                 + x(i,j-1,k) + x(i,j+1,k) 
     $                 + x(i,j,k-1) + x(i,j,k+1) )
                  

                  
               enddo
            enddo
         enddo
      enddo

      return
      end

      subroutine InitArrays(Nk,Nj,Ni, x, y, proc_dim , mycoord, my_offset)
      implicit none

      integer Nk,Nj,Ni
      real*8 x(0:Ni+1,0:Nj+1,0:Nk+1), val
      real*8 y(0:Ni+1,0:Nj+1,0:Nk+1)

      integer proc_dim(1:3), mycoord(1:3), my_offset(1:3)

      integer i,j,k
      integer kstart,kend,jstart,jend,istart,iend

      kstart=0
      jstart=0
      istart=0

      kend=Nk+1
      jend=Nj+1
      iend=Ni+1


      do k = kStart, kEnd
         do j = jStart, jEnd
            do i = iStart, iEnd
               x(i,j,k) = 0.d0
               y(i,j,k) = 0.d0
            enddo
         enddo
      enddo


      if(mycoord(1).eq.0) then
         kstart=0
         jstart=0
         istart=0

         kend=0
         jend=Nj+1
         iend=Ni+1

         val = 10.d0

         call SetArray(Nk,Nj,Ni, x, val,kstart,kend,jstart,jend,istart,iend)

      endif
         
      if(mycoord(2).eq.0) then
         kstart=0
         jstart=0
         istart=0

         kend=Nk+1
         jend=0
         iend=Ni+1

         val = 10.d0

         call SetArray(Nk,Nj,Ni, x, val,kstart,kend,jstart,jend,istart,iend)

      endif
         
      if(mycoord(3).eq.0) then
         kstart=0
         jstart=0
         istart=0

         kend=Nk+1
         jend=Nj+1
         iend=0

         val = 10.d0

         call SetArray(Nk,Nj,Ni, x, val,kstart,kend,jstart,jend,istart,iend)

      endif
         

      if(mycoord(1).eq.proc_dim(1)-1) then
         kstart=Nk+1
         jstart=0
         istart=0

         kend=Nk+1
         jend=Nj+1
         iend=Ni+1

         val = 0.d0

         call SetArray(Nk,Nj,Ni, x, val,kstart,kend,jstart,jend,istart,iend)

      endif
         
      if(mycoord(2).eq.proc_dim(2)-1) then
         kstart=0
         jstart=Nj+1
         istart=0

         kend=Nk+1
         jend=Nj+1
         iend=Ni+1

         val = 0.d0

         call SetArray(Nk,Nj,Ni, x, val,kstart,kend,jstart,jend,istart,iend)

      endif
         
      if(mycoord(3).eq.proc_dim(3)-1) then
         kstart=0
         jstart=0
         istart=Ni+1

         kend=Nk+1
         jend=Nj+1
         iend=Ni+1

         val = 0.d0

         call SetArray(Nk,Nj,Ni, x, val,kstart,kend,jstart,jend,istart,iend)

      endif
         

      return 
      end


      subroutine SetArray(Nk,Nj,Ni, x, val,kstart,kend,jstart,jend,istart,iend)
      implicit none

      integer Nk,Nj,Ni
      real*8 x(0:Ni+1,0:Nj+1,0:Nk+1), val

      integer kstart,kend,jstart,jend,istart,iend
      integer i,j,k

      do k = kStart, kEnd
         do j = jStart, jEnd
            do i = iStart, iEnd
               x(i,j,k) = val
            enddo
         enddo
      enddo
      
      return
      end
