c
c     subroutine to calculate valules for Wilson loops, W, for 'SDLT' shape
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

      subroutine baseup(mx1r,mx1i,mx2r,mx2i,my1r,my1i,my2r,my2i,y1,y2,x1,x2,xdir,ydir,bllr,blli)
c
c     modules
c
      USE L_baryonParam
c
c     no implicit typing
c
      implicit none
c
c     global variables
c
c
c     step links
c
      double precision,dimension(nx,ny,nz,nt,nc,nc)           :: mx1r,mx1i,mx2r,mx2i,my1r,my1i,my2r,my2i
!HPF$ DISTRIBUTE mx1r(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE mx1i(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE mx2r(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE mx2i(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE my1r(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE my1i(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE my2r(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE my2i(*,*,BLOCK,BLOCK,*,*)
      integer                                                 :: y1,y2,x1,x2
      integer                                                 :: xdir,ydir
c
c     bottom links
c
      double precision,dimension(nx,ny,nz,nt,nc,nc)           :: bllr,blli
!HPF$ DISTRIBUTE bllr(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE blli(*,*,BLOCK,BLOCK,*,*)
c
c     local variables
c
      integer                                                 :: ic,jc,kc
      double precision,dimension(nx,ny,nz,nt,nc,nc)           :: lx1r,lx1i,lx2r,lx2i,ly1r,ly1i,ly2r,ly2i
!HPF$ DISTRIBUTE lx1r(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE lx1i(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE lx2r(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE lx2i(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE ly1r(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE ly1i(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE ly2r(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE ly2i(*,*,BLOCK,BLOCK,*,*)
      double precision,dimension(nx,ny,nz,nt,nc,nc)                     :: tmpr,tmpi
!HPF$ DISTRIBUTE tmpr(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE tmpi(*,*,BLOCK,BLOCK,*,*)
c
c     ----------------
c     execution begins
c     ----------------
c
      lx1r = cshift( mx1r, dim=ydir,shift=y1)
      lx1i = cshift( mx1i, dim=ydir,shift=y1)
      lx2r = cshift( cshift( mx2r, dim=ydir,shift=y2), dim=xdir,shift=x1)
      lx2i = cshift( cshift( mx2i, dim=ydir,shift=y2), dim=xdir,shift=x1)
      ly1r = cshift( cshift( my1r, dim=xdir,shift=x1), dim=ydir,shift=y2)
      ly1i = cshift( cshift( my1i, dim=xdir,shift=x1), dim=ydir,shift=y2)
      ly2r = cshift( my2r, dim=xdir,shift=x2)
      ly2i = cshift( my2i, dim=xdir,shift=x2)

      bllr = 0.0d0
      blli = 0.0d0

      do ic = 1,nc
         do jc = 1,nc
            do kc = 1,nc

               bllr(:,:,:,:,ic,jc) = bllr(:,:,:,:,ic,jc)
     &              + lx1r(:,:,:,:,ic,kc) * ly1r(:,:,:,:,jc,kc) + lx1i(:,:,:,:,ic,kc) * ly1i(:,:,:,:,jc,kc)

               blli(:,:,:,:,ic,jc) = blli(:,:,:,:,ic,jc)
     &              - lx1r(:,:,:,:,ic,kc) * ly1i(:,:,:,:,jc,kc) + lx1i(:,:,:,:,ic,kc) * ly1r(:,:,:,:,jc,kc)

            end do
         end do
      end do

      tmpr = 0.0d0
      tmpi = 0.0d0

      do ic = 1,nc
         do jc = 1,nc
            do kc = 1,nc

               tmpr(:,:,:,:,ic,jc) = tmpr(:,:,:,:,ic,jc)
     &              + bllr(:,:,:,:,ic,kc) * lx2r(:,:,:,:,kc,jc) - blli(:,:,:,:,ic,kc) * lx2i(:,:,:,:,kc,jc)

               tmpi(:,:,:,:,ic,jc) = tmpi(:,:,:,:,ic,jc)
     &              + bllr(:,:,:,:,ic,kc) * lx2i(:,:,:,:,kc,jc) + blli(:,:,:,:,ic,kc) * lx2r(:,:,:,:,kc,jc)

            end do
         end do
      end do

      bllr = 0.0d0
      blli = 0.0d0

      do ic = 1,nc
         do jc = 1,nc
            do kc = 1,nc

               bllr(:,:,:,:,ic,jc) = bllr(:,:,:,:,ic,jc)
     &              + tmpr(:,:,:,:,ic,kc) * ly2r(:,:,:,:,jc,kc) + tmpi(:,:,:,:,ic,kc) * ly2i(:,:,:,:,jc,kc)

               blli(:,:,:,:,ic,jc) = blli(:,:,:,:,ic,jc)
     &              - tmpr(:,:,:,:,ic,kc) * ly2i(:,:,:,:,jc,kc) + tmpi(:,:,:,:,ic,kc) * ly2r(:,:,:,:,jc,kc)

            end do
         end do
      end do

      return

      end subroutine baseup


c========================================================================================================

      subroutine basedo(mx1r,mx1i,mx2r,mx2i,my1r,my1i,my2r,my2i,y1,y2,x1,x2,xdir,ydir,brlr,brli)
c
c     modules
c
      USE L_baryonParam
c
c     no implicit typing
c
      implicit none
c
c     global variables
c
c
c     step links
c
      double precision,dimension(nx,ny,nz,nt,nc,nc)           :: mx1r,mx1i,mx2r,mx2i,my1r,my1i,my2r,my2i
!HPF$ DISTRIBUTE mx1r(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE mx1i(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE mx2r(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE mx2i(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE my1r(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE my1i(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE my2r(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE my2i(*,*,BLOCK,BLOCK,*,*)
      integer                                                 :: y1,y2,x1,x2
      integer                                                 :: xdir,ydir
c
c     bottom links
c
      double precision,dimension(nx,ny,nz,nt,nc,nc)           :: brlr,brli
!HPF$ DISTRIBUTE bllr(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE blli(*,*,BLOCK,BLOCK,*,*)
c
c     local variables
c
      integer                                                 :: ic,jc,kc
      double precision,dimension(nx,ny,nz,nt,nc,nc)           :: lx1r,lx1i,lx2r,lx2i,ly1r,ly1i,ly2r,ly2i
!HPF$ DISTRIBUTE lx1r(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE lx1i(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE lx2r(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE lx2i(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE ly1r(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE ly1i(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE ly2r(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE ly2i(*,*,BLOCK,BLOCK,*,*)
      double precision,dimension(nx,ny,nz,nt,nc,nc)                     :: tmpr,tmpi
!HPF$ DISTRIBUTE tmpr(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE tmpi(*,*,BLOCK,BLOCK,*,*)
c
c     ----------------
c     execution begins
c     ----------------
c
      lx1r = cshift( mx1r, dim=ydir,shift=y1)
      lx1i = cshift( mx1i, dim=ydir,shift=y1)
      lx2r = cshift( cshift( mx2r, dim=ydir,shift=y2), dim=xdir,shift=x1)
      lx2i = cshift( cshift( mx2i, dim=ydir,shift=y2), dim=xdir,shift=x1)
      ly1r = cshift( cshift( my1r, dim=xdir,shift=x1), dim=ydir,shift=y1)
      ly1i = cshift( cshift( my1i, dim=xdir,shift=x1), dim=ydir,shift=y1)
      ly2r = cshift( cshift( my2r, dim=xdir,shift=x2), dim=ydir,shift=y2)
      ly2i = cshift( cshift( my2i, dim=xdir,shift=x2), dim=ydir,shift=y2)

      brlr = 0.0d0
      brli = 0.0d0

      do ic = 1,nc
         do jc = 1,nc
            do kc = 1,nc

               brlr(:,:,:,:,ic,jc) = brlr(:,:,:,:,ic,jc)
     &              + lx1r(:,:,:,:,ic,kc) * ly1r(:,:,:,:,kc,jc) - lx1i(:,:,:,:,ic,kc) * ly1i(:,:,:,:,kc,jc)

               brli(:,:,:,:,ic,jc) = brli(:,:,:,:,ic,jc)
     &              + lx1r(:,:,:,:,ic,kc) * ly1i(:,:,:,:,kc,jc) + lx1i(:,:,:,:,ic,kc) * ly1r(:,:,:,:,kc,jc)

            end do
         end do
      end do

      tmpr = 0.0d0
      tmpi = 0.0d0

      do ic = 1,nc
         do jc = 1,nc
            do kc = 1,nc

               tmpr(:,:,:,:,ic,jc) = tmpr(:,:,:,:,ic,jc)
     &              + brlr(:,:,:,:,ic,kc) * lx2r(:,:,:,:,kc,jc) - brli(:,:,:,:,ic,kc) * lx2i(:,:,:,:,kc,jc)

               tmpi(:,:,:,:,ic,jc) = tmpi(:,:,:,:,ic,jc)
     &              + brlr(:,:,:,:,ic,kc) * lx2i(:,:,:,:,kc,jc) + brli(:,:,:,:,ic,kc) * lx2r(:,:,:,:,kc,jc)

            end do
         end do
      end do

      brlr = 0.0d0
      brli = 0.0d0

      do ic = 1,nc
         do jc = 1,nc
            do kc = 1,nc

               brlr(:,:,:,:,ic,jc) = brlr(:,:,:,:,ic,jc)
     &              + tmpr(:,:,:,:,ic,kc) * ly2r(:,:,:,:,kc,jc) - tmpi(:,:,:,:,ic,kc) * ly2i(:,:,:,:,kc,jc)

               brli(:,:,:,:,ic,jc) = brli(:,:,:,:,ic,jc)
     &              + tmpr(:,:,:,:,ic,kc) * ly2i(:,:,:,:,kc,jc) + tmpi(:,:,:,:,ic,kc) * ly2r(:,:,:,:,kc,jc)

            end do
         end do
      end do

      return

      end subroutine basedo

c========================================================================================================

      subroutine delta_loops(xdir,ydir,ur,ui,W)
c
c     modules
c
      USE L_baryonParam
      USE L_epsilonIndex
      USE L_product
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
c     the five cases (1-nYloop) are stored in the second to last place in 'W'
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

      return

      end subroutine delta_loops

      END MODULE L_LOOPS
