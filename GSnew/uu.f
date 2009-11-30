!
!  Collection of routines for link multiplication
!
      MODULE GS_UU

      CONTAINS

      subroutine UU(tor,toi,t1r,t1i,t2r,t2i)

      USE GS_LATTICESIZE

      implicit none
!
      integer,parameter                                           :: nc=3
!
      double precision,dimension(nx,ny,nz,nt,nc,nc),intent(inout) :: tor,toi
!HPF$ DISTRIBUTE tor(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE toi(*,*,BLOCK,BLOCK,*,*)
      double precision,dimension(nx,ny,nz,nt,nc,nc),intent(in)    :: t1r,t1i,t2r,t2i
!HPF$ DISTRIBUTE t1r(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE t1i(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE t2r(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE t2i(*,*,BLOCK,BLOCK,*,*)
!
      integer                                                     :: ic,jc
!
!  Timer Support
!
!      INTEGER start_count, end_count, count_rate
!      REAL    elapsed_time
!
!  Useage
!
!  CALL SYSTEM_CLOCK(start_count, count_rate)
!  CALL SYSTEM_CLOCK(end_count)
!  elapsed_time = REAL(end_count - start_count) / REAL(count_rate)
!
!      CALL SYSTEM_CLOCK(start_count, count_rate)
!
      do ic=1,nc
        do jc=1,nc
          tor(:,:,:,:,ic,jc) = tor(:,:,:,:,ic,jc)
     &      + t1r(:,:,:,:,ic,1) * t2r(:,:,:,:,1,jc) - t1i(:,:,:,:,ic,1) * t2i(:,:,:,:,1,jc)
     &      + t1r(:,:,:,:,ic,2) * t2r(:,:,:,:,2,jc) - t1i(:,:,:,:,ic,2) * t2i(:,:,:,:,2,jc)
     &      + t1r(:,:,:,:,ic,3) * t2r(:,:,:,:,3,jc) - t1i(:,:,:,:,ic,3) * t2i(:,:,:,:,3,jc)
          toi(:,:,:,:,ic,jc) = toi(:,:,:,:,ic,jc)
     &      + t1r(:,:,:,:,ic,1) * t2i(:,:,:,:,1,jc) + t1i(:,:,:,:,ic,1) * t2r(:,:,:,:,1,jc)
     &      + t1r(:,:,:,:,ic,2) * t2i(:,:,:,:,2,jc) + t1i(:,:,:,:,ic,2) * t2r(:,:,:,:,2,jc)
     &      + t1r(:,:,:,:,ic,3) * t2i(:,:,:,:,3,jc) + t1i(:,:,:,:,ic,3) * t2r(:,:,:,:,3,jc)
        end do
      end do
!
!      CALL SYSTEM_CLOCK(end_count)
!      elapsed_time = REAL(end_count - start_count) / REAL(count_rate)
!
!  Excluding assignments
!
!      if (elapsed_time .ne. 0.0) then
!        write(*,'(/,a,f7.3,a)') 'The elapsed time is',elapsed_time,
!     &    ' seconds for one call to UU.'
!        write(*,'(a,f6.3,a)') 'The processing speed is',
!     &    216.0*nx*ny*nz*nt/elapsed_time/1.0e09,' GFLOPS.'
!      end if
!
      end subroutine UU
!
!
!
      subroutine UUdag(tor,toi,t1r,t1i,t2r,t2i)

      USE GS_LATTICESIZE

      implicit none
!
      integer,parameter                                           :: nc=3
!
      double precision,dimension(nx,ny,nz,nt,nc,nc),intent(inout) :: tor,toi
!HPF$ DISTRIBUTE tor(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE toi(*,*,BLOCK,BLOCK,*,*)
      double precision,dimension(nx,ny,nz,nt,nc,nc),intent(in)    :: t1r,t1i,t2r,t2i
!HPF$ DISTRIBUTE t1r(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE t1i(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE t2r(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE t2i(*,*,BLOCK,BLOCK,*,*)
!
      integer                                                     :: ic,jc
!
      do ic=1,nc
        do jc=1,nc
          tor(:,:,:,:,ic,jc) = tor(:,:,:,:,ic,jc)
     &      + t1r(:,:,:,:,ic,1) * t2r(:,:,:,:,jc,1) + t1i(:,:,:,:,ic,1) * t2i(:,:,:,:,jc,1)
     &      + t1r(:,:,:,:,ic,2) * t2r(:,:,:,:,jc,2) + t1i(:,:,:,:,ic,2) * t2i(:,:,:,:,jc,2)
     &      + t1r(:,:,:,:,ic,3) * t2r(:,:,:,:,jc,3) + t1i(:,:,:,:,ic,3) * t2i(:,:,:,:,jc,3)
          toi(:,:,:,:,ic,jc) = toi(:,:,:,:,ic,jc)
     &      - t1r(:,:,:,:,ic,1) * t2i(:,:,:,:,jc,1) + t1i(:,:,:,:,ic,1) * t2r(:,:,:,:,jc,1)
     &      - t1r(:,:,:,:,ic,2) * t2i(:,:,:,:,jc,2) + t1i(:,:,:,:,ic,2) * t2r(:,:,:,:,jc,2)
     &      - t1r(:,:,:,:,ic,3) * t2i(:,:,:,:,jc,3) + t1i(:,:,:,:,ic,3) * t2r(:,:,:,:,jc,3)
        end do
      end do
!
      end subroutine UUdag
!
!
!
      subroutine UdagUdag(tor,toi,t1r,t1i,t2r,t2i)

      USE GS_LATTICESIZE

      implicit none
!
      integer,parameter                                           :: nc=3
!
      double precision,dimension(nx,ny,nz,nt,nc,nc),intent(inout) :: tor,toi
!HPF$ DISTRIBUTE tor(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE toi(*,*,BLOCK,BLOCK,*,*)
      double precision,dimension(nx,ny,nz,nt,nc,nc),intent(in)    :: t1r,t1i,t2r,t2i
!HPF$ DISTRIBUTE t1r(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE t1i(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE t2r(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE t2i(*,*,BLOCK,BLOCK,*,*)
!
      integer                                                     :: ic,jc
!
      do ic=1,nc
        do jc=1,nc
          tor(:,:,:,:,ic,jc) = tor(:,:,:,:,ic,jc)
     &      + t1r(:,:,:,:,1,ic) * t2r(:,:,:,:,jc,1) - t1i(:,:,:,:,1,ic) * t2i(:,:,:,:,jc,1)
     &      + t1r(:,:,:,:,2,ic) * t2r(:,:,:,:,jc,2) - t1i(:,:,:,:,2,ic) * t2i(:,:,:,:,jc,2)
     &      + t1r(:,:,:,:,3,ic) * t2r(:,:,:,:,jc,3) - t1i(:,:,:,:,3,ic) * t2i(:,:,:,:,jc,3)
          toi(:,:,:,:,ic,jc) = toi(:,:,:,:,ic,jc)
     &      - t1r(:,:,:,:,1,ic) * t2i(:,:,:,:,jc,1) - t1i(:,:,:,:,1,ic) * t2r(:,:,:,:,jc,1)
     &      - t1r(:,:,:,:,2,ic) * t2i(:,:,:,:,jc,2) - t1i(:,:,:,:,2,ic) * t2r(:,:,:,:,jc,2)
     &      - t1r(:,:,:,:,3,ic) * t2i(:,:,:,:,jc,3) - t1i(:,:,:,:,3,ic) * t2r(:,:,:,:,jc,3)
        end do
      end do
!
      end subroutine UdagUdag

      END MODULE GS_UU
