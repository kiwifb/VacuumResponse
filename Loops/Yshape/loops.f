c
c     subroutine to calculate valules for Wilson loops, W, for 'Y' shape
c     v.3.11 AK040129
c     v 4.00 FB071129
c
      MODULE L_LOOPS

      CONTAINS

c========================================================================================================

      subroutine y_mirror(xdir,ydir,bxlr,bxli,bllr,blli,brlr,brli,uptr,upti,Res,it,xtop,xbot,yend)
c
c     modules
c
      USE L_baryonParam
      USE L_epsilonIndex
c
c     no implicit typing
c
      implicit none
c
      integer                                                           :: it
      integer                                                           :: xdir,ydir
      integer                                                           :: xtop,xbot,yend
c
c     products of links
c
      double precision,dimension(nx,ny,nz,nt,nc,nc)           		:: uptr,upti
!HPF$ DISTRIBUTE uptr(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE upti(*,*,BLOCK,BLOCK,*,*)
c
c     bottom and top links (built up iteratively)
c
      double precision,dimension(nx,ny,nz,nt,nc,nc)           		:: bxlr,bxli
!HPF$ DISTRIBUTE bxlr(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE bxli(*,*,BLOCK,BLOCK,*,*)
      double precision,dimension(nx,ny,nz,nt,nc,nc)           		:: bllr,blli
!HPF$ DISTRIBUTE bllr(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE blli(*,*,BLOCK,BLOCK,*,*)
      double precision,dimension(nx,ny,nz,nt,nc,nc)           		:: brlr,brli
!HPF$ DISTRIBUTE brlr(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE brli(*,*,BLOCK,BLOCK,*,*)
      double precision,dimension(nx,ny,nz,nt,nc,nc)           		:: tlr,tli
!HPF$ DISTRIBUTE tlr(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE tli(*,*,BLOCK,BLOCK,*,*)
c
c     Mirror bottom links
c
      double precision,dimension(nx,ny,nz,nt,nc,nc)           		:: mbxlr,mbxli
!HPF$ DISTRIBUTE mbxlr(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE mbxli(*,*,BLOCK,BLOCK,*,*)
      double precision,dimension(nx,ny,nz,nt,nc,nc)           		:: mbllr,mblli
!HPF$ DISTRIBUTE mbllr(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE mblli(*,*,BLOCK,BLOCK,*,*)
      double precision,dimension(nx,ny,nz,nt,nc,nc)           		:: mbrlr,mbrli
!HPF$ DISTRIBUTE mbrlr(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE mbrli(*,*,BLOCK,BLOCK,*,*)
c
c     shifted time links
c
      double precision,dimension(nx,ny,nz,nt,nc,nc)           		:: ustr,usti
!HPF$ DISTRIBUTE ustr(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE usti(*,*,BLOCK,BLOCK,*,*)
c
c     staples/temporary arrays
c
      double precision,dimension(nx,ny,nz,nt,nc,nc)           		:: sxlr,sxli
!HPF$ DISTRIBUTE sxlr(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE sxli(*,*,BLOCK,BLOCK,*,*)
      double precision,dimension(nx,ny,nz,nt,nc,nc)           		:: sllr,slli
!HPF$ DISTRIBUTE sllr(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE slli(*,*,BLOCK,BLOCK,*,*)
      double precision,dimension(nx,ny,nz,nt,nc,nc)           		:: srlr,srli
!HPF$ DISTRIBUTE srlr(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE srli(*,*,BLOCK,BLOCK,*,*)
c
c     temporary arrays
c
      double precision,dimension(nx,ny,nz,nt,nc,nc)           		:: tmpr,tmpi
!HPF$ DISTRIBUTE tmpr(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE tmpi(*,*,BLOCK,BLOCK,*,*)
      double precision,dimension(nx,ny,nz,nt)                           :: Res
!HPF$ DISTRIBUTE Res(*,*,BLOCK,BLOCK)
c
c     local variables
c
      integer                                                 		:: ic,jc,kc
c
c     Mirror top (-xtop,0)
c
      mbxlr = cshift(bxlr,dim=xdir,shift=-xtop)
      mbxli = cshift(bxli,dim=xdir,shift=-xtop)

      tmpr = 0.0d0
      tmpi = 0.0d0
      ustr = cshift(uptr,dim=xdir,shift=-xtop)
      usti = cshift(upti,dim=xdir,shift=-xtop)

      do ic = 1,nc
         do jc = 1,nc
            do kc = 1,nc

               tmpr(:,:,:,:,ic,jc) = tmpr(:,:,:,:,ic,jc)
     &              + mbxlr(:,:,:,:,kc,ic) * ustr(:,:,:,:,kc,jc) + mbxli(:,:,:,:,kc,ic) * usti(:,:,:,:,kc,jc)

               tmpi(:,:,:,:,ic,jc) = tmpi(:,:,:,:,ic,jc)
     &              + mbxlr(:,:,:,:,kc,ic) * usti(:,:,:,:,kc,jc) - mbxli(:,:,:,:,kc,ic) * ustr(:,:,:,:,kc,jc)

            end do
         end do
      end do

      sxlr = 0.0d0
      sxli = 0.0d0
      tlr = cshift(mbxlr,dim=that,shift=it)
      tli = cshift(mbxli,dim=that,shift=it)

      do ic = 1,nc
         do jc = 1,nc
            do kc = 1,nc

               sxlr(:,:,:,:,ic,jc) = sxlr(:,:,:,:,ic,jc)
     &              + tmpr(:,:,:,:,ic,kc) * tlr(:,:,:,:,kc,jc) - tmpi(:,:,:,:,ic,kc) * tli(:,:,:,:,kc,jc)

               sxli(:,:,:,:,ic,jc) = sxli(:,:,:,:,ic,jc)
     &              + tmpr(:,:,:,:,ic,kc) * tli(:,:,:,:,kc,jc) + tmpi(:,:,:,:,ic,kc) * tlr(:,:,:,:,kc,jc)

            end do
         end do
      end do
c
c     Mirror left, i.e., quark at (-xbot,-yend)
c
      mbllr = cshift(cshift(bllr,dim=xdir,shift=-xbot),dim=ydir,shift=-yend)
      mblli = cshift(cshift(blli,dim=xdir,shift=-xbot),dim=ydir,shift=-yend)

      tmpr = 0.0d0
      tmpi = 0.0d0
      ustr = cshift(cshift(uptr,dim=xdir,shift=-xbot),dim=ydir,shift=-yend)
      usti = cshift(cshift(upti,dim=xdir,shift=-xbot),dim=ydir,shift=-yend)

      do ic = 1,nc
         do jc = 1,nc
            do kc = 1,nc

               tmpr(:,:,:,:,ic,jc) = tmpr(:,:,:,:,ic,jc)
     &              + mbllr(:,:,:,:,kc,ic) * ustr(:,:,:,:,kc,jc) + mblli(:,:,:,:,kc,ic) * usti(:,:,:,:,kc,jc)

               tmpi(:,:,:,:,ic,jc) = tmpi(:,:,:,:,ic,jc)
     &              + mbllr(:,:,:,:,kc,ic) * usti(:,:,:,:,kc,jc) - mblli(:,:,:,:,kc,ic) * ustr(:,:,:,:,kc,jc)

            end do
         end do
      end do

      sllr = 0.0d0
      slli = 0.0d0
      tlr = cshift(mbllr,dim=that,shift=it)
      tli = cshift(mblli,dim=that,shift=it)

      do ic = 1,nc
         do jc = 1,nc
            do kc = 1,nc

               sllr(:,:,:,:,ic,jc) = sllr(:,:,:,:,ic,jc)
     &              + tmpr(:,:,:,:,ic,kc) * tlr(:,:,:,:,kc,jc) - tmpi(:,:,:,:,ic,kc) * tli(:,:,:,:,kc,jc)

               slli(:,:,:,:,ic,jc) = slli(:,:,:,:,ic,jc)
     &              + tmpr(:,:,:,:,ic,kc) * tli(:,:,:,:,kc,jc) + tmpi(:,:,:,:,ic,kc) * tlr(:,:,:,:,kc,jc)

            end do
         end do
      end do
c
c     Mirror right, i.e., quark at (-xbot,+yend)
c
      mbrlr = cshift(cshift(brlr,dim=xdir,shift=-xbot),dim=ydir,shift= yend)
      mbrli = cshift(cshift(brli,dim=xdir,shift=-xbot),dim=ydir,shift= yend)

      tmpr = 0.0d0
      tmpi = 0.0d0
      ustr = cshift(cshift(uptr,dim=xdir,shift=-xbot),dim=ydir,shift= yend)
      usti = cshift(cshift(upti,dim=xdir,shift=-xbot),dim=ydir,shift= yend)

      do ic = 1,nc
         do jc = 1,nc
            do kc = 1,nc

               tmpr(:,:,:,:,ic,jc) = tmpr(:,:,:,:,ic,jc)
     &              + mbrlr(:,:,:,:,kc,ic) * ustr(:,:,:,:,kc,jc) + mbrli(:,:,:,:,kc,ic) * usti(:,:,:,:,kc,jc)

               tmpi(:,:,:,:,ic,jc) = tmpi(:,:,:,:,ic,jc)
     &              + mbrlr(:,:,:,:,kc,ic) * usti(:,:,:,:,kc,jc) - mbrli(:,:,:,:,kc,ic) * ustr(:,:,:,:,kc,jc)

            end do
         end do
      end do

      srlr = 0.0d0
      srli = 0.0d0
      tlr = cshift(mbrlr,dim=that,shift=it)
      tli = cshift(mbrli,dim=that,shift=it)

      do ic = 1,nc
         do jc = 1,nc
            do kc = 1,nc

               srlr(:,:,:,:,ic,jc) = srlr(:,:,:,:,ic,jc)
     &              + tmpr(:,:,:,:,ic,kc) * tlr(:,:,:,:,kc,jc) - tmpi(:,:,:,:,ic,kc) * tli(:,:,:,:,kc,jc)

               srli(:,:,:,:,ic,jc) = srli(:,:,:,:,ic,jc)
     &              + tmpr(:,:,:,:,ic,kc) * tli(:,:,:,:,kc,jc) + tmpi(:,:,:,:,ic,kc) * tlr(:,:,:,:,kc,jc)

            end do
         end do
      end do

      include 'res.f'

      return

      end subroutine y_mirror

c========================================================================================================

      subroutine y_loops(xdir,ydir,ur,ui,W)
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
      integer                                                                   :: xdir,ydir
      double precision,dimension(nx,ny,nz,nt,mu,nc,nc)            	 	:: ur,ui
!HPF$ DISTRIBUTE ur(*,*,BLOCK,BLOCK,*,*,*)
!HPF$ DISTRIBUTE ui(*,*,BLOCK,BLOCK,*,*,*)
c
c     the five cases (0-nYloop) are stored in the second to last place in 'W'
c
      double precision,dimension(nx,ny,nz,nt,0:nloop,nL(that))          :: W
!HPF$ DISTRIBUTE W(*,*,BLOCK,BLOCK,*,*)
c
      include 'loop_declarations.f'
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
      call boxm1p2(xdir,ydir,ur,ui,up1x2r,up1x2i)
      call boxm1m2(xdir,ydir,ur,ui,do1x2r,do1x2i)
      call boxm2p3(xdir,ydir,ur,ui,up2x3r,up2x3i)
      call boxm2m3(xdir,ydir,ur,ui,do2x3r,do2x3i)
c
c     ---------------------------------
c     start of Wilson loop calculations
c     ---------------------------------
c
c     it do loop
c
      do it = 1, nL(that)

         call product(ur,ui,that,uptr,upti,it)
c
c     CASE 0 - three quarks at (0,0)
c
         include 'Czero.f'

         W(:,:,:,:,0,it) = Res(:,:,:,:) / 6.0d0
c
c     CASE 1 - quarks at ( 1,0), (-1,-1) and (-1, 1)
c
         include 'Cone.f'

         W(:,:,:,:,1,it) = Res(:,:,:,:) / 6.0d0
c
c     CASE 2 - quarks at ( 2,0), (-1,-2) and (-1, 2)
c
         include 'Ctwo.f'

         W(:,:,:,:,2,it) = Res(:,:,:,:) / 6.0d0
c
c     CASE 3 - quarks at ( 3,0), (-1,-2) and (-1, 2)
c
         include 'Cthree.f'

         W(:,:,:,:,3,it) = Res(:,:,:,:) / 6.0d0
c
c     CASE 4 - quarks at ( 3,0), (-2,-3) and (-2, 3)
c
         include 'Cfour.f'

         W(:,:,:,:,4,it) = Res(:,:,:,:) / 6.0d0
c
c     CASE 5 - quarks at ( 4,0), (-3,-4) and (-3, 4)
c
         include 'Cfive.f'

         W(:,:,:,:,5,it) = Res(:,:,:,:) / 6.0d0
c
c     CASE 6 - quarks at ( 5,0), (-4,-5) and (-4, 5)
c
         include 'Csix.f'

         W(:,:,:,:,6,it) = Res(:,:,:,:) / 6.0d0
c
c     CASE 7 - quarks at ( 7,0), (-4,-6) and (-4, 6)
c
         include 'Cseven.f'

         W(:,:,:,:,7,it) = Res(:,:,:,:) / 6.0d0

#ifdef __S16__
c
c     CASE 8 - quarks at ( 8,0), (-4,-7) and (-4, 7)
c
         include 'Ceight.f'

         W(:,:,:,:,8,it) = Res(:,:,:,:) / 6.0d0
c
c     CASE 9 - quarks at ( 9,0), (-5,-8) and (-5, 8)
c
         include 'Cnine.f'

         W(:,:,:,:,9,it) = Res(:,:,:,:) / 6.0d0
#endif

      end do

      return

      end subroutine y_loops

c========================================================================================================

      subroutine y_loops_mirror(xdir,ydir,ur,ui,Wp,Wm)
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
      double precision,dimension(nx,ny,nz,nt,mu,nc,nc)            	 	:: ur,ui
!HPF$ DISTRIBUTE ur(*,*,BLOCK,BLOCK,*,*,*)
!HPF$ DISTRIBUTE ui(*,*,BLOCK,BLOCK,*,*,*)
c
c     the five cases (0-nYloop) are stored in the second to last place in 'W'
c
      double precision,dimension(nx,ny,nz,nt,0:nloop,nL(that))         :: Wp,Wm
!HPF$ DISTRIBUTE Wp(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE Wm(*,*,BLOCK,BLOCK,*,*)
c
      include 'loop_declarations.f'
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
      call boxm1p2(xdir,ydir,ur,ui,up1x2r,up1x2i)
      call boxm1m2(xdir,ydir,ur,ui,do1x2r,do1x2i)
      call boxm2p3(xdir,ydir,ur,ui,up2x3r,up2x3i)
      call boxm2m3(xdir,ydir,ur,ui,do2x3r,do2x3i)
c
c     ---------------------------------
c     start of Wilson loop calculations
c     ---------------------------------
c
c     it do loop
c
      do it = 1, nL(that)

         call product(ur,ui,that,uptr,upti,it)
c
c     CASE 0 - three quarks at (0,0)
c
         include 'Czero.f'

         Wp(:,:,:,:,0,it) = Res(:,:,:,:) / 6.0d0
         Wm(:,:,:,:,0,it) = Res(:,:,:,:) / 6.0d0
c
c     CASE 1 - quarks at ( 1,0), (-1,-1) and (-1, 1)
c     Mirror - quarks at (-1,0), ( 1, 1) and ( 1,-1)
c
         include 'Cone.f'

         Wp(:,:,:,:,1,it) = Res(:,:,:,:) / 6.0d0

         call y_mirror(xdir,ydir,bxlr,bxli,bllr,blli,brlr,brli,uptr,upti,Res,it,1,-1,1)

         Wm(:,:,:,:,1,it) = Res(:,:,:,:) / 6.0d0
c
c     CASE 2 - quarks at ( 2,0), (-1,-2) and (-1, 2)
c     Mirror - quarks at (-2,0), ( 1, 2) and ( 1,-2)
c
         include 'Ctwo.f'

         Wp(:,:,:,:,2,it) = Res(:,:,:,:) / 6.0d0

         call y_mirror(xdir,ydir,bxlr,bxli,bllr,blli,brlr,brli,uptr,upti,Res,it,2,-1,2)

         Wm(:,:,:,:,2,it) = Res(:,:,:,:) / 6.0d0
c
c     CASE 3 - quarks at ( 3,0), (-1,-2) and (-1, 2)
c     Mirror - quarks at (-3,0), ( 1, 2) and ( 1,-2)
c
         include 'Cthree.f'

         Wp(:,:,:,:,3,it) = Res(:,:,:,:) / 6.0d0
c
         call y_mirror(xdir,ydir,bxlr,bxli,bllr,blli,brlr,brli,uptr,upti,Res,it,3,-1,2)

         Wm(:,:,:,:,3,it) = Res(:,:,:,:) / 6.0d0
c
c     CASE 4 - quarks at ( 3,0), (-2,-3) and (-2, 3)
c     Mirror - quarks at (-3,0), ( 2, 3) and ( 3,-3)
c
         include 'Cfour.f'

         Wp(:,:,:,:,4,it) = Res(:,:,:,:) / 6.0d0
c
         call y_mirror(xdir,ydir,bxlr,bxli,bllr,blli,brlr,brli,uptr,upti,Res,it,3,-2,3)

         Wm(:,:,:,:,4,it) = Res(:,:,:,:) / 6.0d0
c
c     CASE 5 - quarks at ( 4,0), (-3,-4) and (-3, 4)
c     Mirror - quarks at (-4,0), ( 3, 4) and ( 3,-4)
c
         include 'Cfive.f'

         Wp(:,:,:,:,5,it) = Res(:,:,:,:) / 6.0d0
c
         call y_mirror(xdir,ydir,bxlr,bxli,bllr,blli,brlr,brli,uptr,upti,Res,it,4,-3,4)

         Wm(:,:,:,:,5,it) = Res(:,:,:,:) / 6.0d0
c
c     CASE 6 - quarks at ( 5,0), (-4,-5) and (-4, 5)
c     Mirror - quarks at (-5,0), ( 4, 5) and ( 4,-5)
c
         include 'Csix.f'

         Wp(:,:,:,:,6,it) = Res(:,:,:,:) / 6.0d0
c
         call y_mirror(xdir,ydir,bxlr,bxli,bllr,blli,brlr,brli,uptr,upti,Res,it,5,-4,5)

         Wm(:,:,:,:,6,it) = Res(:,:,:,:) / 6.0d0
c
c     CASE 7 - quarks at ( 7,0), (-4,-6) and (-4, 6)
c     Mirror - quarks at (-7,0), ( 4, 6) and ( 4,-6)
c
         include 'Cseven.f'

         Wp(:,:,:,:,7,it) = Res(:,:,:,:) / 6.0d0
c
         call y_mirror(xdir,ydir,bxlr,bxli,bllr,blli,brlr,brli,uptr,upti,Res,it,7,-4,6)

         Wm(:,:,:,:,7,it) = Res(:,:,:,:) / 6.0d0
#ifdef __S16__
c
c     CASE 8 - quarks at ( 8,0), (-4,-7) and (-4, 7)
c     Mirror - quarks at (-8,0), ( 4, 7) and ( 4,-7)
c
         include 'Ceight.f'

         Wp(:,:,:,:,8,it) = Res(:,:,:,:) / 6.0d0
c
         call y_mirror(xdir,ydir,bxlr,bxli,bllr,blli,brlr,brli,uptr,upti,Res,it,8,-4,7)

         Wm(:,:,:,:,8,it) = Res(:,:,:,:) / 6.0d0
c
c     CASE 9 - quarks at ( 9,0), (-5,-8) and (-5, 8)
c     Mirror - quarks at (-9,0), ( 5, 8) and ( 5,-8)
c
         include 'Cnine.f'

         Wp(:,:,:,:,9,it) = Res(:,:,:,:) / 6.0d0
c
         call y_mirror(xdir,ydir,bxlr,bxli,bllr,blli,brlr,brli,uptr,upti,Res,it,9,-5,8)

         Wm(:,:,:,:,9,it) = Res(:,:,:,:) / 6.0d0
#endif

      end do

      return

      end subroutine y_loops_mirror

      END MODULE L_LOOPS
