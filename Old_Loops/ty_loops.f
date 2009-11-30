c
c     subroutines to calculate valules for Wilson loops, W, for 'T' shape
c     corresponding to 'Y' shape loops.
c     v.3.11 AK040129
c
      MODULE L_TYLOOPS

      include'Yfiles/yLoopSize.h'

      CONTAINS

c========================================================================================================

      subroutine ty_mirror(xdir,ydir,bxlr,bxli,sllr,slli,srlr,srli,uptr,upti,Res,it,xtop)
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
      integer                                                           :: xtop
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
      double precision,dimension(nx,ny,nz,nt,nc,nc)           		:: tlr,tli
!HPF$ DISTRIBUTE tlr(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE tli(*,*,BLOCK,BLOCK,*,*)
c
c     Mirror bottom links 
c
      double precision,dimension(nx,ny,nz,nt,nc,nc)           		:: mbxlr,mbxli
!HPF$ DISTRIBUTE mbxlr(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE mbxli(*,*,BLOCK,BLOCK,*,*)
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
c     No need for mirror left or right we are studying T shapes!
c     We already have the staples
c
      include 'TYfiles/res.f'

      return

      end subroutine ty_mirror

c========================================================================================================

      subroutine ty_loops_mirror(xdir,ydir,ur,ui,Wp,Wm)
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
      integer                                                           :: xdir,ydir
      double precision,dimension(nx,ny,nz,nt,mu,nc,nc)         	 	:: ur,ui
!HPF$ DISTRIBUTE ur(*,*,BLOCK,BLOCK,*,*,*)
!HPF$ DISTRIBUTE ui(*,*,BLOCK,BLOCK,*,*,*)
c
c     the five cases (0-nYloop) are stored in the second to last place in 'W'
c
      double precision,dimension(nx,ny,nz,nt,0:nYloop,nL(that))         :: Wp,Wm
!HPF$ DISTRIBUTE Wp(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE Wm(*,*,BLOCK,BLOCK,*,*)
c
      include 'TYfiles/loop_declarations.f'
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
         include 'TYfiles/Czero.f'

         Wp(:,:,:,:,0,it) = Res(:,:,:,:) / 6.0d0
         Wm(:,:,:,:,0,it) = Res(:,:,:,:) / 6.0d0
c
c     CASE 1 - quarks at ( 1,0), (-1,-1) and (-1, 1) = ( 2,0), (0, 1) and (0,-1)
c     Mirror - quarks at (-1,0), ( 1, 1) and ( 1,-1) = (-2,0), (0, 1) and (0,-1)
c
         include 'TYfiles/Cone.f'

         Wp(:,:,:,:,1,it) = Res(:,:,:,:) / 6.0d0

         call ty_mirror(xdir,ydir,bxlr,bxli,sllr,slli,srlr,srli,uptr,upti,Res,it,2)

         Wm(:,:,:,:,1,it) = Res(:,:,:,:) / 6.0d0
c
c     CASE 2 - quarks at ( 2,0), (-1,-2) and (-1, 2) = ( 3,0), (0, 2) and (0,-2)
c     Mirror - quarks at (-2,0), ( 1, 2) and ( 1,-2) = (-3,0), (0, 2) and (0,-2)
c
         include 'TYfiles/Ctwo.f'

         Wp(:,:,:,:,2,it) = Res(:,:,:,:) / 6.0d0

         call ty_mirror(xdir,ydir,bxlr,bxli,sllr,slli,srlr,srli,uptr,upti,Res,it,3)

         Wm(:,:,:,:,2,it) = Res(:,:,:,:) / 6.0d0
c
c     CASE 3 - quarks at ( 3,0), (-1,-2) and (-1, 2) = ( 4,0), (0, 2) and (0,-2)
c     Mirror - quarks at (-3,0), ( 1, 2) and ( 1,-2) = (-4,0), (0, 2) and (0,-2)
c
         include 'TYfiles/Cthree.f'

         Wp(:,:,:,:,3,it) = Res(:,:,:,:) / 6.0d0
c
         call ty_mirror(xdir,ydir,bxlr,bxli,sllr,slli,srlr,srli,uptr,upti,Res,it,4)

         Wm(:,:,:,:,3,it) = Res(:,:,:,:) / 6.0d0
c
c     CASE 4 - quarks at ( 3,0), (-2,-3) and (-2, 3) = ( 5,0), (0, 3) and (0,-3)
c     Mirror - quarks at (-3,0), ( 2, 3) and ( 3,-3) = (-5,0), (0, 3) and (0,-3)
c
         include 'TYfiles/Cfour.f'

         Wp(:,:,:,:,4,it) = Res(:,:,:,:) / 6.0d0
c
         call ty_mirror(xdir,ydir,bxlr,bxli,sllr,slli,srlr,srli,uptr,upti,Res,it,5)

         Wm(:,:,:,:,4,it) = Res(:,:,:,:) / 6.0d0
c
c     CASE 5 - quarks at ( 4,0), (-3,-4) and (-3, 4) = ( 7,0), (0, 4) and (0,-4)
c     Mirror - quarks at (-4,0), ( 3, 4) and ( 3,-4) = (-7,0), (0, 4) and (0,-4)
c
         include 'TYfiles/Cfive.f'

         Wp(:,:,:,:,5,it) = Res(:,:,:,:) / 6.0d0
c
         call ty_mirror(xdir,ydir,bxlr,bxli,sllr,slli,srlr,srli,uptr,upti,Res,it,7)

         Wm(:,:,:,:,5,it) = Res(:,:,:,:) / 6.0d0
c
c     CASE 6 - quarks at ( 5,0), (-4,-5) and (-4, 5) = ( 9,0), (0, 5) and (0,-5)
c     Mirror - quarks at (-5,0), ( 4, 5) and ( 4,-5) = (-9,0), (0, 5) and (0,-5)
c
         include 'TYfiles/Csix.f'

         Wp(:,:,:,:,6,it) = Res(:,:,:,:) / 6.0d0
c
         call ty_mirror(xdir,ydir,bxlr,bxli,sllr,slli,srlr,srli,uptr,upti,Res,it,9)

         Wm(:,:,:,:,6,it) = Res(:,:,:,:) / 6.0d0
c
c     CASE 7 - quarks at ( 7,0), (-4,-6) and (-4, 6) = ( 11,0), (0, 6) and (0,-6)
c     Mirror - quarks at (-7,0), ( 4, 6) and ( 4,-6) = (-11,0), (0, 6) and (0,-6)
c
         include 'TYfiles/Cseven.f'

         Wp(:,:,:,:,7,it) = Res(:,:,:,:) / 6.0d0
c
         call ty_mirror(xdir,ydir,bxlr,bxli,sllr,slli,srlr,srli,uptr,upti,Res,it,11)

         Wm(:,:,:,:,7,it) = Res(:,:,:,:) / 6.0d0

         if(nYloop == 9) then
c
c     CASE 8 - quarks at ( 8,0), (-4,-7) and (-4, 7) = ( 12,0), (0, 7) and (0,-7)
c     Mirror - quarks at (-8,0), ( 4, 7) and ( 4,-7) = (-12,0), (0, 7) and (0,-7)
c
            include 'TYfiles/Ceight.f'

            Wp(:,:,:,:,8,it) = Res(:,:,:,:) / 6.0d0
c
            call ty_mirror(xdir,ydir,bxlr,bxli,sllr,slli,srlr,srli,uptr,upti,Res,it,12)

            Wm(:,:,:,:,8,it) = Res(:,:,:,:) / 6.0d0

c
c     CASE 9 - quarks at ( 9,0), (-5,-8) and (-5, 8) = ( 14,0), (0, 8) and (0,-8)
c     Mirror - quarks at (-9,0), ( 5, 8) and ( 5,-8) = (-14,0), (0, 8) and (0,-8)
c
            include 'TYfiles/Cnine.f'

            Wp(:,:,:,:,9,it) = Res(:,:,:,:) / 6.0d0
c
            call ty_mirror(xdir,ydir,bxlr,bxli,sllr,slli,srlr,srli,uptr,upti,Res,it,14)

            Wm(:,:,:,:,9,it) = Res(:,:,:,:) / 6.0d0

         endif

      end do

      return

      end subroutine ty_loops_mirror

      END MODULE L_TYLOOPS
