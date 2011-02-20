c
c     subroutine to calculate valules for Wilson loops, W, for 'VY' shape
!        adapted from 'Delta' shape
c     1.0 20060110 FB
!     2.0 FB071129 - ported to use boxes.f and hopefully correct a bug causing gauge
!        invariance problem in shape 7.
c
      MODULE L_LOOPS

      CONTAINS

c========================================================================================================

      subroutine link_conjug(ur,ui)
c
c     modules
c
      USE GS_LATTICESIZE
      USE L_baryonParam
c
c     no implicit typing
c
      implicit none
c
c     global variables
c
      double precision,dimension(nx,ny,nz,nt,nc,nc)                  :: ur,ui
!HPF$ DISTRIBUTE ur(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE ui(*,*,BLOCK,BLOCK,*,*)
c
c     local variables
c
      double precision,dimension(nx,ny,nz,nt,nc,nc)                     :: tmpr,tmpi
!HPF$ DISTRIBUTE tmpr(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE tmpi(*,*,BLOCK,BLOCK,*,*)
      integer                                                           :: ic,jc

      do ic = 1, nc
         do jc = 1, nc

            tmpr(:,:,:,:,ic,jc) = ur(:,:,:,:,jc,ic)
            tmpi(:,:,:,:,ic,jc) =-ui(:,:,:,:,jc,ic)

         end do
      end do

      ur = tmpr
      ui = tmpi

      return

      end subroutine link_conjug

c========================================================================================================

      subroutine delta_loops(xdir,ydir,ur,ui,W)
c
c     modules
c
      USE L_baryonParam
      USE L_epsilonIndex
      USE L_product
      USE L_BOXES
c
c     no implicit typing
c
      implicit none
c
c     global variables
c
#include"loopsize.f"
      integer                                                           :: xdir,ydir
      double precision,dimension(nx,ny,nz,nt,mu,nc,nc)                  :: ur,ui
!HPF$ DISTRIBUTE ur(*,*,BLOCK,BLOCK,*,*,*)
!HPF$ DISTRIBUTE ui(*,*,BLOCK,BLOCK,*,*,*)
c
c     intermediate links
c
      double precision,dimension(nx,ny,nz,nt,nc,nc)           :: ilr,ili
!HPF$ DISTRIBUTE ilr(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE ili(*,*,BLOCK,BLOCK,*,*)
      double precision,dimension(nx,ny,nz,nt,nc,nc)           :: irr,iri
!HPF$ DISTRIBUTE irr(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE iri(*,*,BLOCK,BLOCK,*,*)
      double precision,dimension(nx,ny,nz,nt,nc,nc)           :: ixr, ixi
!HPF$ DISTRIBUTE ixr(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE ixi(*,*,BLOCK,BLOCK,*,*)
      double precision,dimension(nx,ny,nz,nt,nc,nc)           :: uyr, uyi, dyr, dyi
!HPF$ DISTRIBUTE uyr(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE uyi(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE dyr(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE dyi(*,*,BLOCK,BLOCK,*,*)
c
c     all the cases (1-nYloop) are stored in the second to last place in 'W'
c
      double precision,dimension(nx,ny,nz,nt,0:nloop,nL(that))          :: W
!HPF$ DISTRIBUTE W(*,*,BLOCK,BLOCK,*,*)
      integer                                                    :: time1, time2, count_rate
      integer                                                    :: timea, timeb
c
      include'loop_declarations.f'
c
c     ----------------
c     execution begins
c     ----------------
c
c     -----------------------
c     make the diagonal links
c     -----------------------
c
      call boxm1p1(xdir,ydir,ur,ui,up1x1r,up1x1i)
      call boxm1m1(xdir,ydir,ur,ui,do1x1r,do1x1i)
      call boxm2p1(xdir,ydir,ur,ui,up2x1r,up2x1i)
      call boxm2m1(xdir,ydir,ur,ui,do2x1r,do2x1i)
      call boxm3p2(xdir,ydir,ur,ui,up3x2r,up3x2i)
      call boxm3m2(xdir,ydir,ur,ui,do3x2r,do3x2i)
c
c     --------------------------
c     prepare the straight bits
c     --------------------------
c
      ixr(:,:,:,:,:,:) = ur(:,:,:,:,xdir,:,:)
      ixi(:,:,:,:,:,:) = ui(:,:,:,:,xdir,:,:)
      call link_conjug(ixr,ixi)
      ixr = cshift(ixr, dim=xdir,shift=-1)
      ixi = cshift(ixi, dim=xdir,shift=-1)

      uyr(:,:,:,:,:,:) = ur(:,:,:,:,ydir,:,:)
      uyi(:,:,:,:,:,:) = ui(:,:,:,:,ydir,:,:)

      dyr = cshift(uyr, dim=ydir,shift=-1)
      dyi = cshift(uyi, dim=ydir,shift=-1)
      call link_conjug(dyr,dyi)
c
c     ---------------------------------
c     start of Wilson loop calculations
c     ---------------------------------
c
c     compute the time link loop
c
      do it = 1, nL(that)

         call product(ur,ui,that,tmp1r,tmp1i,it)
         tlr(:,:,:,:,:,:,it) = tmp1r(:,:,:,:,:,:)
         tli(:,:,:,:,:,:,it) = tmp1i(:,:,:,:,:,:)

      end do
c
c     CASE 1 - quarks at ( 1,0), (-1,-1) and (-1, 1)
c
      call system_clock(time1,count_rate)
      include'Cone.f'
      call system_clock(time2)
      write(*,'(a,f15.8,a)') 'loops 1: ',(real(time2-time1)/real(count_rate)),' s'
c
c     CASE 2 - quarks at ( 2,0), (-1,-2) and (-1, 2)
c
      include'Ctwo.f'
      call system_clock(time1)
      write(*,'(a,f15.8,a)') 'loops 2: ',(real(time1-time2)/real(count_rate)),' s'
c
c     CASE 3 - quarks at ( 3,0), (-1,-2) and (-1, 2)
c
      include'Cthree.f'
      call system_clock(time2)
      write(*,'(a,f15.8,a)') 'loops 3: ',(real(time2-time1)/real(count_rate)),' s'
c
c     CASE 4 - quarks at ( 3,0), (-2,-3) and (-2, 3)
c
      include'Cfour.f'
      call system_clock(time1)
      write(*,'(a,f15.8,a)') 'loops 4: ',(real(time1-time2)/real(count_rate)),' s'
c
c     CASE 5 - quarks at ( 4,0), (-3,-4) and (-3, 4)
c
      include'Cfive.f'
      call system_clock(time2)
      write(*,'(a,f15.8,a)') 'loops 5: ',(real(time2-time1)/real(count_rate)),' s'
c
c     CASE 6 - quarks at ( 5,0), (-4,-5) and (-4, 5)
c
      include'Csix.f'
      call system_clock(time1)
      write(*,'(a,f15.8,a)') 'loops 6: ',(real(time1-time2)/real(count_rate)),' s'
c
c     CASE 7 - quarks at ( 7,0), (-4,-6) and (-4, 6)
c
      include'Cseven.f'
      call system_clock(time2)
      write(*,'(a,f15.8,a)') 'loops 7: ',(real(time2-time1)/real(count_rate)),' s'
c      if(nloop == 9) then
c
c     CASE 8 - quarks at ( 8,0), (-4,-7) and (-4, 7)
c
c         include'Ceight.f'
c
c         W(:,:,:,:,8,it) = Res(:,:,:,:) / 216.0d0
c
c     CASE 9 - quarks at ( 9,0), (-5,-8) and (-5, 8)
c
c         include'Cnine.f'
c
c         W(:,:,:,:,9,it) = Res(:,:,:,:) / 216.0d0
c
c      endif

      return

      end subroutine delta_loops

      END MODULE L_LOOPS
