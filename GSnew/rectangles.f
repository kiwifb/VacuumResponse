!
!     cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!     computes the product of links for the action associated with a link in the
!     xhat direction. This computes the product of the six rectangular plaquettes
!     associated with the improved action.
!
!     Authors: Frederi! D.R. Bonnet & Derek B. Leinweber
!     Last update: Derek B. Leinweber, February 2001
!
      MODULE GS_RECTANGLES

      CONTAINS

      subroutine rectangles(ur,ui,rectr,recti,xhat,local)

      USE GS_LATTICESIZE
      USE GS_UU

      implicit none

!     global variable

      integer,parameter                                       :: nc=3               !sigma,color
      integer,parameter                                       :: mu=4               !direction

      integer                                                 :: xhat
      logical                                                 :: local

      double precision,dimension(nx,ny,nz,nt,nc,nc)           :: rectr,recti        !rectangular staples
!HPF$ DISTRIBUTE rectr(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE recti(*,*,BLOCK,BLOCK,*,*)
      double precision,dimension(nx,ny,nz,nt,mu,nc,nc)        :: ur,ui
!HPF$ DISTRIBUTE ur(*,*,BLOCK,BLOCK,*,*,*)
!HPF$ DISTRIBUTE ui(*,*,BLOCK,BLOCK,*,*,*)

!     local variables, temporary product variables

      integer,dimension(3)                                    :: yhat               !direction vector
!HPF$ DISTRIBUTE yhat(*)
      double precision,dimension(nx,ny,nz,nt,nc,nc)           :: t1r,t1i,t2r,t2i    !return variables
!HPF$ DISTRIBUTE t1r(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE t1i(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE t2r(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE t2i(*,*,BLOCK,BLOCK,*,*)
      double precision,dimension(nx,ny,nz,nt,nc,nc)           :: usr,usi
!HPF$ DISTRIBUTE usr(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE usi(*,*,BLOCK,BLOCK,*,*)

      integer                                                 :: imu,i             !counters
!c
!!  Timer Support
!c
!      INTEGER start_count, end_count, count_rate
!      REAL elapsed_time
!
!  Useage
!
!  CALL SYSTEM_CLOCK(start_count, count_rate)
!  CALL SYSTEM_CLOCK(end_count)
!  elapsed_time = REAL(end_count - start_count) / REAL(count_rate)

!     starting of the execution commands

!     setting up the yhat array
!     the yhat(1)=yhat,yhat(2)=zhat and yhat(3)=that when xhat eq 1
!     the yhat(1)=zhat,yhat(2)=that and yhat(3)=xhat when xhat eq 2
!     the yhat(1)=that,yhat(2)=xhat and yhat(3)=yhat when xhat eq 3
!     the yhat(1)=xhat,yhat(2)=yhat and yhat(3)=zhat when xhat eq 4

      do imu=1,mu-1
         yhat(imu) = 0
      end do

      do imu=1,mu-1
         yhat(imu)    = mod( xhat + imu ,4 )
         if(yhat(imu) .eq. 0) then
            yhat(imu) = 4
         end if
      end do

!      CALL SYSTEM_CLOCK(start_count, count_rate)

!     calculation of the link products in the positive plaquette.
!     starting point in the xhat direction with the maskrect

!     first calculate the positive forward 2a x a rectangle

      rectr = 0.0d0
      recti = 0.0d0

      do i=1,mu-1                     !loop over the yhat array
         t1r = 0.0d0
         t1i = 0.0d0
         t2r = cshift(ur(:,:,:,:,xhat,:,:),dim=xhat,shift=1)
         t2i = cshift(ui(:,:,:,:,xhat,:,:),dim=xhat,shift=1)
         usr = cshift(ur(:,:,:,:,yhat(i),:,:),dim=xhat,shift=2)
         usi = cshift(ui(:,:,:,:,yhat(i),:,:),dim=xhat,shift=2)
         call UU(t1r,t1i,t2r,t2i,usr,usi)

         t2r = 0.0d0
         t2i = 0.0d0
         usr = cshift(cshift(ur(:,:,:,:,xhat,:,:),dim=yhat(i),shift=1),dim=xhat,shift=1)
         usi = cshift(cshift(ui(:,:,:,:,xhat,:,:),dim=yhat(i),shift=1),dim=xhat,shift=1)
         call UUdag(t2r,t2i,t1r,t1i,usr,usi)

         t1r = 0.0d0
         t1i = 0.0d0
         usr = cshift(ur(:,:,:,:,xhat,:,:),dim=yhat(i),shift=1)
         usi = cshift(ui(:,:,:,:,xhat,:,:),dim=yhat(i),shift=1)
         call UUdag(t1r,t1i,t2r,t2i,usr,usi)

         usr = ur(:,:,:,:,yhat(i),:,:)
         usi = ui(:,:,:,:,yhat(i),:,:)
         call UUdag(rectr,recti,t1r,t1i,usr,usi)

!     now calculating the positive upper rectangles a x 2a

         t1r = 0.0d0
         t1i = 0.0d0
         t2r = cshift(ur(:,:,:,:,yhat(i),:,:),dim=xhat,shift=1)
         t2i = cshift(ui(:,:,:,:,yhat(i),:,:),dim=xhat,shift=1)
         usr = cshift(cshift(ur(:,:,:,:,yhat(i),:,:),dim=xhat,shift=1),dim=yhat(i),shift=1)
         usi = cshift(cshift(ui(:,:,:,:,yhat(i),:,:),dim=xhat,shift=1),dim=yhat(i),shift=1)
         call UU(t1r,t1i,t2r,t2i,usr,usi)

         t2r = 0.0d0
         t2i = 0.0d0
         usr = cshift(ur(:,:,:,:,xhat,:,:),dim=yhat(i),shift=2)
         usi = cshift(ui(:,:,:,:,xhat,:,:),dim=yhat(i),shift=2)
         call UUdag(t2r,t2i,t1r,t1i,usr,usi)

         t1r = 0.0d0
         t1i = 0.0d0
         usr = cshift(ur(:,:,:,:,yhat(i),:,:),dim=yhat(i),shift=1)
         usi = cshift(ui(:,:,:,:,yhat(i),:,:),dim=yhat(i),shift=1)
         call UUdag(t1r,t1i,t2r,t2i,usr,usi)

         usr = ur(:,:,:,:,yhat(i),:,:)
         usi = ui(:,:,:,:,yhat(i),:,:)
         call UUdag(rectr,recti,t1r,t1i,usr,usi)

!     .true. when we calculate the full staple
!     .false. when we calculte impbar where
!     only the upper 3 plaquette are needed

         if(local) then

!     now calculate the positive backward 2a x a retangle

            t1r = 0.0d0
            t1i = 0.0d0
            t2r = cshift(ur(:,:,:,:,yhat(i),:,:),dim=xhat,shift=1)
            t2i = cshift(ui(:,:,:,:,yhat(i),:,:),dim=xhat,shift=1)
            usr = cshift(ur(:,:,:,:,xhat,:,:),dim=yhat(i),shift=1)
            usi = cshift(ui(:,:,:,:,xhat,:,:),dim=yhat(i),shift=1)
            call UUdag(t1r,t1i,t2r,t2i,usr,usi)

            t2r = 0.0d0
            t2i = 0.0d0
            usr = cshift(cshift(ur(:,:,:,:,xhat,:,:),dim=yhat(i),shift=1),dim=xhat,shift=-1)
            usi = cshift(cshift(ui(:,:,:,:,xhat,:,:),dim=yhat(i),shift=1),dim=xhat,shift=-1)
            call UUdag(t2r,t2i,t1r,t1i,usr,usi)

            t1r = 0.0d0
            t1i = 0.0d0
            usr = cshift(ur(:,:,:,:,yhat(i),:,:),dim=xhat,shift=-1)
            usi = cshift(ui(:,:,:,:,yhat(i),:,:),dim=xhat,shift=-1)
            call UUdag(t1r,t1i,t2r,t2i,usr,usi)

            usr = cshift(ur(:,:,:,:,xhat,:,:),dim=xhat,shift=-1)
            usi = cshift(ui(:,:,:,:,xhat,:,:),dim=xhat,shift=-1)
            call UU(rectr,recti,t1r,t1i,usr,usi)

!     calculation of the link products in the negative plaquette
!     starting point in the x direction

!     the negative xy,xz and xt contour when xhat eq 1
!     the negative yz,yt and yx contour when xhat eq 2
!     the negative zt,zx and zy contour when xhat eq 3
!     the negative tx,ty and tz contour when xhat eq 4


!     first calculate the negative forward 2a x a rectangle

            t1r = 0.0d0
            t1i = 0.0d0
            t2r = cshift(ur(:,:,:,:,xhat,:,:),dim=xhat,shift=1)
            t2i = cshift(ui(:,:,:,:,xhat,:,:),dim=xhat,shift=1)
            usr = cshift(cshift(ur(:,:,:,:,yhat(i),:,:),dim=xhat,shift=2),dim=yhat(i),shift=-1)
            usi = cshift(cshift(ui(:,:,:,:,yhat(i),:,:),dim=xhat,shift=2),dim=yhat(i),shift=-1)
            call UUdag(t1r,t1i,t2r,t2i,usr,usi)

            t2r = 0.0d0
            t2i = 0.0d0
            usr = cshift(cshift(ur(:,:,:,:,xhat,:,:),dim=yhat(i),shift=-1),dim=xhat,shift=1)
            usi = cshift(cshift(ui(:,:,:,:,xhat,:,:),dim=yhat(i),shift=-1),dim=xhat,shift=1)
            call UUdag(t2r,t2i,t1r,t1i,usr,usi)

            t1r = 0.0d0
            t1i = 0.0d0
            usr = cshift(ur(:,:,:,:,xhat,:,:),dim=yhat(i),shift=-1)
            usi = cshift(ui(:,:,:,:,xhat,:,:),dim=yhat(i),shift=-1)
            call UUdag(t1r,t1i,t2r,t2i,usr,usi)

            usr = cshift(ur(:,:,:,:,yhat(i),:,:),dim=yhat(i),shift=-1)
            usi = cshift(ui(:,:,:,:,yhat(i),:,:),dim=yhat(i),shift=-1)
            call UU(rectr,recti,t1r,t1i,usr,usi)

!     calculate the negative backward 2a x a rectangle

            t1r = 0.0d0
            t1i = 0.0d0
            t2r = cshift(cshift(ur(:,:,:,:,yhat(i),:,:),dim=yhat(i),shift=-1),dim=xhat,shift=1)
            t2i = cshift(cshift(ui(:,:,:,:,yhat(i),:,:),dim=yhat(i),shift=-1),dim=xhat,shift=1)
            usr = cshift(ur(:,:,:,:,xhat,:,:),dim=yhat(i),shift=-1)
            usi = cshift(ui(:,:,:,:,xhat,:,:),dim=yhat(i),shift=-1)
            call UdagUdag(t1r,t1i,t2r,t2i,usr,usi)

            t2r = 0.0d0
            t2i = 0.0d0
            usr = cshift(cshift(ur(:,:,:,:,xhat,:,:),dim=yhat(i),shift=-1),dim=xhat,shift=-1)
            usi = cshift(cshift(ui(:,:,:,:,xhat,:,:),dim=yhat(i),shift=-1),dim=xhat,shift=-1)
            call UUdag(t2r,t2i,t1r,t1i,usr,usi)

            t1r = 0.0d0
            t1i = 0.0d0
            usr = cshift(cshift(ur(:,:,:,:,yhat(i),:,:),dim=xhat,shift=-1),dim=yhat(i),shift=-1)
            usi = cshift(cshift(ui(:,:,:,:,yhat(i),:,:),dim=xhat,shift=-1),dim=yhat(i),shift=-1)
            call UU(t1r,t1i,t2r,t2i,usr,usi)

            usr = cshift(ur(:,:,:,:,xhat,:,:),dim=xhat,shift=-1)
            usi = cshift(ui(:,:,:,:,xhat,:,:),dim=xhat,shift=-1)
            call UU(rectr,recti,t1r,t1i,usr,usi)

!     now calculating the negative lower rectangles a x 2a

            t1r = 0.0d0
            t1i = 0.0d0
            t2r = cshift(cshift(ur(:,:,:,:,yhat(i),:,:),dim=xhat,shift=1),dim=yhat(i),shift=-1)
            t2i = cshift(cshift(ui(:,:,:,:,yhat(i),:,:),dim=xhat,shift=1),dim=yhat(i),shift=-1)
            usr = cshift(cshift(ur(:,:,:,:,yhat(i),:,:),dim=xhat,shift=1),dim=yhat(i),shift=-2)
            usi = cshift(cshift(ui(:,:,:,:,yhat(i),:,:),dim=xhat,shift=1),dim=yhat(i),shift=-2)
            call UdagUdag(t1r,t1i,t2r,t2i,usr,usi)

            t2r = 0.0d0
            t2i = 0.0d0
            usr = cshift(ur(:,:,:,:,xhat,:,:),dim=yhat(i),shift=-2)
            usi = cshift(ui(:,:,:,:,xhat,:,:),dim=yhat(i),shift=-2)
            call UUdag(t2r,t2i,t1r,t1i,usr,usi)

            t1r = 0.0d0
            t1i = 0.0d0
            usr = cshift(ur(:,:,:,:,yhat(i),:,:),dim=yhat(i),shift=-2)
            usi = cshift(ui(:,:,:,:,yhat(i),:,:),dim=yhat(i),shift=-2)
            call UU(t1r,t1i,t2r,t2i,usr,usi)

            usr = cshift(ur(:,:,:,:,yhat(i),:,:),dim=yhat(i),shift=-1)
            usi = cshift(ui(:,:,:,:,yhat(i),:,:),dim=yhat(i),shift=-1)
            call UU(rectr,recti,t1r,t1i,usr,usi)

         end if
      end do

!      CALL SYSTEM_CLOCK(end_count)
!      elapsed_time = REAL(end_count - start_count) / REAL(count_rate)
!c
!!  Excluding assignments
!c
!      if (local) then
!        write(*,'(/,a,f7.3,a)') 'The elapsed time is',elapsed_time,
!     &    ' seconds for 6 rectangles.'
!        write(*,'(a,f6.3,a)') 'The processing speed is',
!     &    6.0*2592.0*nx*ny*nz*nt/elapsed_time/1.0e09,' GFLOPS.'
!      else
!        write(*,'(/,a,f7.3,a)') 'The elapsed time is',elapsed_time,
!     &    ' seconds for 2 rectangles.'
!        write(*,'(a,f6.3,a)') 'The processing speed is',
!     &    2.0*2592.0*nx*ny*nz*nt/elapsed_time/1.0e09,' GFLOPS.'
!      end if

      return

      end subroutine rectangles

      END MODULE GS_RECTANGLES
