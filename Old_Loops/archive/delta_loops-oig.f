c
c     subroutine to calculate valules for Wilson loops, W, for 'Delta' shape
c     1.0 20060110 FB
c
      MODULE L_DELTALOOPS

      include'DELTAfiles/DeltaLoopSize.h'

      CONTAINS

c
c     Create boxes of links for the Delta-shape loops
c     They are oriented differently from the Y-shape loops
c     Here is the space orientation:
c
c        y ^
c          |
c          o--> x
c
c     box1x1 creates the two following boxes:
c
c     boxp       boxm
c
c     o->o       o->o
c     |  |  and  ^  ^
c     V  V       |  |
c     o->o       o->o
c
c     Note that they are not transposed of each other
c
c     box2x1 creates the two following boxes:
c
c      boxp          boxm
c
c     o->o->o       o->o->o
c     |     |  and  ^     ^
c     V     V       |     |
c     o->o->o       o->o->o
c

c========================================================================================================

      subroutine box1x1(xdir,ydir,ur,ui,boxpr,boxpi,boxmr,boxmi)

      USE L_baryonParam

      IMPLICIT NONE

      double precision,dimension(nx,ny,nz,nt,mu,nc,nc)        :: ur,ui !links
!HPF$ DISTRIBUTE ur(*,*,BLOCK,BLOCK,*,*,*)
!HPF$ DISTRIBUTE ui(*,*,BLOCK,BLOCK,*,*,*)
      double precision,dimension(nx,ny,nz,nt,nc,nc)           :: boxpr,boxpi,boxmr,boxmi
!HPF$ DISTRIBUTE boxpr(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE boxpi(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE boxmr(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE boxmi(*,*,BLOCK,BLOCK,*,*)

      double precision,dimension(nx,ny,nz,nt,nc,nc)           :: ux1r,ux1i,ux2r,ux2i,uy1r,uy1i,uy2r,uy2i
!HPF$ DISTRIBUTE ux1r(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE ux1i(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE ux2r(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE ux2i(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE uy1r(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE uy1i(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE uy2r(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE uy2i(*,*,BLOCK,BLOCK,*,*)

      integer                                                 :: xdir,ydir
      integer                                                 :: ic,jc,kc      

      boxpr = 0.d0
      boxpi = 0.d0
      boxmr = 0.d0
      boxmi = 0.d0

      ux1r(:,:,:,:,:,:) = ur(:,:,:,:,xdir,:,:)
      ux1i(:,:,:,:,:,:) = ui(:,:,:,:,xdir,:,:)

      uy1r(:,:,:,:,:,:) = cshift(ur(:,:,:,:,ydir,:,:), dim = ydir, shift =-1)
      uy1i(:,:,:,:,:,:) = cshift(ui(:,:,:,:,ydir,:,:), dim = ydir, shift =-1)

      ux2r(:,:,:,:,:,:) = cshift( ux1r(:,:,:,:,:,:), dim = ydir, shift =-1)
      ux2i(:,:,:,:,:,:) = cshift( ux1i(:,:,:,:,:,:), dim = ydir, shift =-1)

      uy2r(:,:,:,:,:,:) = cshift( uy1r(:,:,:,:,:,:), dim = xdir, shift = 1)
      uy2i(:,:,:,:,:,:) = cshift( uy1i(:,:,:,:,:,:), dim = xdir, shift = 1)

      do ic = 1, nc
         do jc = 1, nc
            do kc = 1, nc

               boxpr(:,:,:,:,ic,jc) = boxpr(:,:,:,:,ic,jc) + 
     &              ux1r(:,:,:,:,ic,kc)*uy2r(:,:,:,:,jc,kc) + ux1i(:,:,:,:,ic,kc)*uy2i(:,:,:,:,jc,kc) +
     &              uy1r(:,:,:,:,kc,ic)*ux2r(:,:,:,:,kc,jc) + uy1i(:,:,:,:,kc,ic)*ux2i(:,:,:,:,kc,jc)

               boxpi(:,:,:,:,ic,jc) = boxpi(:,:,:,:,ic,jc) - 
     &              ux1r(:,:,:,:,ic,kc)*uy2i(:,:,:,:,jc,kc) + ux1i(:,:,:,:,ic,kc)*uy2r(:,:,:,:,jc,kc) +
     &              uy1r(:,:,:,:,kc,ic)*ux2i(:,:,:,:,kc,jc) - uy1i(:,:,:,:,kc,ic)*ux2r(:,:,:,:,kc,jc)

            end do
         end do
      end do

      uy1r(:,:,:,:,:,:) = ur(:,:,:,:,ydir,:,:)
      uy1i(:,:,:,:,:,:) = ui(:,:,:,:,ydir,:,:)

      ux2r(:,:,:,:,:,:) = cshift( ux1r(:,:,:,:,:,:), dim = ydir, shift = 1)
      ux2i(:,:,:,:,:,:) = cshift( ux1i(:,:,:,:,:,:), dim = ydir, shift = 1)

      uy2r(:,:,:,:,:,:) = cshift( uy1r(:,:,:,:,:,:), dim = xdir, shift = 1)
      uy2i(:,:,:,:,:,:) = cshift( uy1i(:,:,:,:,:,:), dim = xdir, shift = 1)

      do ic = 1, nc
         do jc = 1, nc
            do kc = 1, nc

               boxmr(:,:,:,:,ic,jc) = boxmr(:,:,:,:,ic,jc) + 
     &              ux1r(:,:,:,:,ic,kc)*uy2r(:,:,:,:,kc,jc) - ux1i(:,:,:,:,ic,kc)*uy2i(:,:,:,:,kc,jc) +
     &              uy1r(:,:,:,:,ic,kc)*ux2r(:,:,:,:,kc,jc) - uy1i(:,:,:,:,ic,kc)*ux2i(:,:,:,:,kc,jc)

               boxmi(:,:,:,:,ic,jc) = boxmi(:,:,:,:,ic,jc) + 
     &              ux1r(:,:,:,:,ic,kc)*uy2i(:,:,:,:,kc,jc) + ux1i(:,:,:,:,ic,kc)*uy2r(:,:,:,:,kc,jc) +
     &              uy1r(:,:,:,:,ic,kc)*ux2i(:,:,:,:,kc,jc) + uy1i(:,:,:,:,ic,kc)*ux2r(:,:,:,:,kc,jc)

            end do
         end do
      end do

      return

      end subroutine box1x1

c========================================================================================================

      subroutine box2x1(xdir,ydir,ur,ui,boxpr,boxpi,boxmr,boxmi)

      USE L_baryonParam

      IMPLICIT NONE

      double precision,dimension(nx,ny,nz,nt,mu,nc,nc)        :: ur,ui !links
!HPF$ DISTRIBUTE ur(*,*,BLOCK,BLOCK,*,*,*)
!HPF$ DISTRIBUTE ui(*,*,BLOCK,BLOCK,*,*,*)
      double precision,dimension(nx,ny,nz,nt,nc,nc)           :: boxpr,boxpi,boxmr,boxmi
!HPF$ DISTRIBUTE boxpr(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE boxpi(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE boxmr(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE boxmi(*,*,BLOCK,BLOCK,*,*)

      double precision,dimension(nx,ny,nz,nt,nc,nc)           :: ux1r,ux1i,ux2r,ux2i,uy1r,uy1i,uy2r,uy2i
!HPF$ DISTRIBUTE ux1r(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE ux1i(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE ux2r(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE ux2i(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE uy1r(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE uy1i(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE uy2r(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE uy2i(*,*,BLOCK,BLOCK,*,*)

      integer                                                 :: xdir,ydir
      integer                                                 :: ic,jc,kc

      boxpr = 0.d0
      boxpi = 0.d0
      boxmr = 0.d0
      boxmi = 0.d0
      ux1r  = 0.d0
      ux1i  = 0.d0

      ux2r(:,:,:,:,:,:) = cshift( ur(:,:,:,:,xdir,:,:), dim = xdir, shift = 1)
      ux2i(:,:,:,:,:,:) = cshift( ui(:,:,:,:,xdir,:,:), dim = xdir, shift = 1)

      do ic = 1, nc
         do jc = 1, nc
            do kc= 1, nc

               ux1r(:,:,:,:,ic,jc) = ux1r(:,:,:,:,ic,jc) + 
     &              ur(:,:,:,:,xdir,ic,kc)*ux2r(:,:,:,:,kc,jc) - ui(:,:,:,:,xdir,ic,kc)*ux2i(:,:,:,:,kc,jc) 

               ux1i(:,:,:,:,ic,jc) = ux1i(:,:,:,:,ic,jc) + 
     &              ur(:,:,:,:,xdir,ic,kc)*ux2i(:,:,:,:,kc,jc) + ui(:,:,:,:,xdir,ic,kc)*ux2r(:,:,:,:,kc,jc) 

            end do
         end do
      end do

      uy1r(:,:,:,:,:,:) = cshift( ur(:,:,:,:,ydir,:,:), dim = ydir, shift =-1)
      uy1i(:,:,:,:,:,:) = cshift( ui(:,:,:,:,ydir,:,:), dim = ydir, shift =-1)

      ux2r(:,:,:,:,:,:) = cshift( ux1r(:,:,:,:,:,:), dim = ydir, shift =-1)
      ux2i(:,:,:,:,:,:) = cshift( ux1i(:,:,:,:,:,:), dim = ydir, shift =-1)

      uy2r(:,:,:,:,:,:) = cshift( uy1r(:,:,:,:,:,:), dim = xdir, shift = 2)
      uy2i(:,:,:,:,:,:) = cshift( uy1i(:,:,:,:,:,:), dim = xdir, shift = 2)

      do ic = 1, nc
         do jc = 1, nc
            do kc = 1, nc

               boxpr(:,:,:,:,ic,jc) = boxpr(:,:,:,:,ic,jc) + 
     &              ux1r(:,:,:,:,ic,kc)*uy2r(:,:,:,:,jc,kc) + ux1i(:,:,:,:,ic,kc)*uy2i(:,:,:,:,jc,kc) +
     &              uy1r(:,:,:,:,kc,ic)*ux2r(:,:,:,:,kc,jc) + uy1i(:,:,:,:,kc,ic)*ux2i(:,:,:,:,kc,jc)

               boxpi(:,:,:,:,ic,jc) = boxpi(:,:,:,:,ic,jc) - 
     &              ux1r(:,:,:,:,ic,kc)*uy2i(:,:,:,:,jc,kc) + ux1i(:,:,:,:,ic,kc)*uy2r(:,:,:,:,jc,kc) +
     &              uy1r(:,:,:,:,kc,ic)*ux2i(:,:,:,:,kc,jc) - uy1i(:,:,:,:,kc,ic)*ux2r(:,:,:,:,kc,jc)

            end do
         end do
      end do

      uy1r(:,:,:,:,:,:) = ur(:,:,:,:,ydir,:,:)
      uy1i(:,:,:,:,:,:) = ui(:,:,:,:,ydir,:,:)

      ux2r(:,:,:,:,:,:) = cshift( ux1r(:,:,:,:,:,:), dim = ydir, shift = 1)
      ux2i(:,:,:,:,:,:) = cshift( ux1i(:,:,:,:,:,:), dim = ydir, shift = 1)

      uy2r(:,:,:,:,:,:) = cshift( uy1r(:,:,:,:,:,:), dim = xdir, shift = 2)
      uy2i(:,:,:,:,:,:) = cshift( uy1i(:,:,:,:,:,:), dim = xdir, shift = 2)

      do ic = 1, nc
         do jc = 1, nc
            do kc = 1, nc

               boxmr(:,:,:,:,ic,jc) = boxmr(:,:,:,:,ic,jc) + 
     &              ux1r(:,:,:,:,ic,kc)*uy2r(:,:,:,:,kc,jc) - ux1i(:,:,:,:,ic,kc)*uy2i(:,:,:,:,kc,jc) +
     &              uy1r(:,:,:,:,ic,kc)*ux2r(:,:,:,:,kc,jc) - uy1i(:,:,:,:,ic,kc)*ux2i(:,:,:,:,kc,jc)

               boxmi(:,:,:,:,ic,jc) = boxmi(:,:,:,:,ic,jc) + 
     &              ux1r(:,:,:,:,ic,kc)*uy2i(:,:,:,:,kc,jc) + ux1i(:,:,:,:,ic,kc)*uy2r(:,:,:,:,kc,jc) +
     &              uy1r(:,:,:,:,ic,kc)*ux2i(:,:,:,:,kc,jc) + uy1i(:,:,:,:,ic,kc)*ux2r(:,:,:,:,kc,jc)
            end do
         end do
      end do

      return

      end subroutine box2x1

c========================================================================================================

      subroutine box2x3(xdir,ydir,bp1r,bp1i,bm1r,bm1i,bp2r,bp2i,bm2r,bm2i,boxpr,boxpi,boxmr,boxmi)

      USE L_baryonParam

      IMPLICIT NONE

      double precision,dimension(nx,ny,nz,nt,nc,nc)           :: bp1r,bp1i,bm1r,bm1i
!HPF$ DISTRIBUTE bp1r(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE bp1i(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE bm1r(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE bm1i(*,*,BLOCK,BLOCK,*,*)
      double precision,dimension(nx,ny,nz,nt,nc,nc)           :: bp2r,bp2i,bm2r,bm2i
!HPF$ DISTRIBUTE bp2r(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE bp2i(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE bm2r(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE bm2i(*,*,BLOCK,BLOCK,*,*)
      double precision,dimension(nx,ny,nz,nt,nc,nc)           :: boxpr,boxpi,boxmr,boxmi
!HPF$ DISTRIBUTE boxpr(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE boxpi(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE boxmr(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE boxmi(*,*,BLOCK,BLOCK,*,*)

      double precision,dimension(nx,ny,nz,nt,nc,nc)           :: sb1r,sb1i,sb2r,sb2i
!HPF$ DISTRIBUTE sb1r(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE sb1i(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE sb2r(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE sb2i(*,*,BLOCK,BLOCK,*,*)

      integer                                                 :: xdir,ydir
      integer                                                 :: ic,jc,kc

      boxpr=0.d0
      boxpi=0.d0

      sb1r = cshift(cshift(bp1r,dim=xdir,shift= 2),dim=ydir,shift=-1)
      sb1i = cshift(cshift(bp1i,dim=xdir,shift= 2),dim=ydir,shift=-1)

      sb2r = cshift(cshift(bp2r,dim=xdir,shift= 1),dim=ydir,shift=-1)
      sb2i = cshift(cshift(bp2i,dim=xdir,shift= 1),dim=ydir,shift=-1)

      do ic = 1, nc
         do jc = 1, nc
            do kc = 1, nc

               boxpr(:,:,:,:,ic,jc) = boxpr(:,:,:,:,ic,jc) + 
     &              bp1r(:,:,:,:,ic,kc)*sb2r(:,:,:,:,kc,jc) - bp1i(:,:,:,:,ic,kc)*sb2i(:,:,:,:,kc,jc) +
     &              bp2r(:,:,:,:,ic,kc)*sb1r(:,:,:,:,kc,jc) - bp2i(:,:,:,:,ic,kc)*sb1i(:,:,:,:,kc,jc)

               boxpi(:,:,:,:,ic,jc) = boxpi(:,:,:,:,ic,jc) + 
     &              bp1r(:,:,:,:,ic,kc)*sb2i(:,:,:,:,kc,jc) + bp1i(:,:,:,:,ic,kc)*sb2r(:,:,:,:,kc,jc) +
     &              bp2r(:,:,:,:,ic,kc)*sb1i(:,:,:,:,kc,jc) + bp2i(:,:,:,:,ic,kc)*sb1r(:,:,:,:,kc,jc)

            end do
         end do
      end do

      boxmr=0.d0
      boxmi=0.d0

      sb1r = cshift(cshift(bm1r,dim=xdir,shift= 2),dim=ydir,shift= 1)
      sb1i = cshift(cshift(bm1i,dim=xdir,shift= 2),dim=ydir,shift= 1)

      sb2r = cshift(cshift(bm2r,dim=xdir,shift= 1),dim=ydir,shift= 1)
      sb2i = cshift(cshift(bm2i,dim=xdir,shift= 1),dim=ydir,shift= 1)

      do ic = 1, nc
         do jc = 1, nc
            do kc = 1, nc

               boxmr(:,:,:,:,ic,jc) = boxmr(:,:,:,:,ic,jc) + 
     &              bm1r(:,:,:,:,ic,kc)*sb2r(:,:,:,:,kc,jc) - bm1i(:,:,:,:,ic,kc)*sb2i(:,:,:,:,kc,jc) +
     &              bm2r(:,:,:,:,ic,kc)*sb1r(:,:,:,:,kc,jc) - bm2i(:,:,:,:,ic,kc)*sb1i(:,:,:,:,kc,jc)

               boxmi(:,:,:,:,ic,jc) = boxmi(:,:,:,:,ic,jc) + 
     &              bm1r(:,:,:,:,ic,kc)*sb2i(:,:,:,:,kc,jc) + bm1i(:,:,:,:,ic,kc)*sb2r(:,:,:,:,kc,jc) +
     &              bm2r(:,:,:,:,ic,kc)*sb1i(:,:,:,:,kc,jc) + bm2i(:,:,:,:,ic,kc)*sb1r(:,:,:,:,kc,jc)

            end do
         end do
      end do

      return

      end subroutine box2x3

c========================================================================================================

      subroutine delta_mirror(xdir,ydir,bxlr,bxli,mbllr,mblli,mbrlr,mbrli,tlr,tli,tl1r,tl1i,Res,it,xtop,ytop)
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
      integer                                                           :: xtop,ytop
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
c
c     time links
c
      double precision,dimension(nx,ny,nz,nt,nc,nc)                    :: tlr,tli
!HPF$ DISTRIBUTE tlr(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE tli(*,*,BLOCK,BLOCK,*,*)
      double precision,dimension(nx,ny,nz,nt,nc,nc)                    :: tl1r,tl1i,tl2r,tl2i
!HPF$ DISTRIBUTE tl1r(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE tl1i(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE tl2r(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE tl2i(*,*,BLOCK,BLOCK,*,*)
c
c     Mirror diagonal links 
c
      double precision,dimension(nx,ny,nz,nt,nc,nc)           		:: mbllr,mblli
!HPF$ DISTRIBUTE mbllr(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE mblli(*,*,BLOCK,BLOCK,*,*)
      double precision,dimension(nx,ny,nz,nt,nc,nc)           		:: mbrlr,mbrli
!HPF$ DISTRIBUTE mbrlr(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE mbrli(*,*,BLOCK,BLOCK,*,*)
c
c     top links (time shifted bottom links)
c
      double precision,dimension(nx,ny,nz,nt,nc,nc)                     :: sbxlr,sbxli
!HPF$ DISTRIBUTE sbxlr(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE sbxli(*,*,BLOCK,BLOCK,*,*)
      double precision,dimension(nx,ny,nz,nt,nc,nc)                     :: sbllr,sblli
!HPF$ DISTRIBUTE sbllr(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE sblli(*,*,BLOCK,BLOCK,*,*)
      double precision,dimension(nx,ny,nz,nt,nc,nc)                     :: sbrlr,sbrli
!HPF$ DISTRIBUTE sbrlr(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE sbrli(*,*,BLOCK,BLOCK,*,*)
c
c     temporary arrays
c
      double precision,dimension(nx,ny,nz,nt)                          :: tmp1r,tmp1i,tmp2r,tmp2i,tmp3r,tmp3i
!HPF$ DISTRIBUTE tmp1r(*,*,BLOCK,BLOCK)
!HPF$ DISTRIBUTE tmp1i(*,*,BLOCK,BLOCK)
!HPF$ DISTRIBUTE tmp2r(*,*,BLOCK,BLOCK)
!HPF$ DISTRIBUTE tmp2i(*,*,BLOCK,BLOCK)
!HPF$ DISTRIBUTE tmp3r(*,*,BLOCK,BLOCK)
!HPF$ DISTRIBUTE tmp3i(*,*,BLOCK,BLOCK)
      double precision,dimension(nx,ny,nz,nt)                           :: Res
!HPF$ DISTRIBUTE Res(*,*,BLOCK,BLOCK)
c
c     local variables
c
      integer                                                 		:: ic,jc,kc
c
c     base is unchanged
c
c     Mirror top diagonal
c
      bllr = cshift(cshift(mbrlr,dim=xdir,shift=-xtop),dim=ydir,shift= ytop)
      blli = cshift(cshift(mbrli,dim=xdir,shift=-xtop),dim=ydir,shift= ytop)
c
c     Mirror bottom diagonal
c
      brlr = cshift(cshift(mbllr,dim=xdir,shift=-xtop),dim=ydir,shift=-ytop)
      brli = cshift(cshift(mblli,dim=xdir,shift=-xtop),dim=ydir,shift=-ytop)
c
c     Mirror time shifted links
c
      tl2r = cshift( cshift( tlr, dim=ydir,shift= ytop),dim=xdir,shift=-xtop)
      tl2i = cshift( cshift( tli, dim=ydir,shift= ytop),dim=xdir,shift=-xtop)

      include 'DELTAfiles/res.f'

      return

      end subroutine delta_mirror

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
      integer                                                           :: xdir,ydir
      double precision,dimension(nx,ny,nz,nt,mu,nc,nc)         	 	:: ur,ui
!HPF$ DISTRIBUTE ur(*,*,BLOCK,BLOCK,*,*,*)
!HPF$ DISTRIBUTE ui(*,*,BLOCK,BLOCK,*,*,*)
c
c     the five cases (1-nYloop) are stored in the second to last place in 'W'
c
      double precision,dimension(nx,ny,nz,nt,0:nYloop,nL(that))         :: W
!HPF$ DISTRIBUTE W(*,*,BLOCK,BLOCK,*,*)
c
      include 'DELTAfiles/loop_declarations.f'
c
c     ----------------
c     execution begins
c     ----------------
c
c     -----------------------
c     make the diagonal links
c     -----------------------
c
      call box1x1(xdir,ydir,ur,ui,ld1r,ld1i,rd1r,rd1i)
      call box2x1(xdir,ydir,ur,ui,ld2r,ld2i,rd2r,rd2i)
      call box2x3(xdir,ydir,ld1r,ld1i,rd1r,rd1i,ld2r,ld2i,rd2r,rd2i,ld3r,ld3i,rd3r,rd3i)
c
c     ---------------------------------
c     start of Wilson loop calculations
c     ---------------------------------
c
c     it do loop
c
      do it = 1, nL(that)

         call product(ur,ui,that,tlr,tli,it)
c
c     CASE 1 - quarks at ( 1,0), (-1,-1) and (-1, 1)
c
         include 'DELTAfiles/Cone.f'

         W(:,:,:,:,1,it) = Res(:,:,:,:) / 216.0d0
c
c     CASE 2 - quarks at ( 2,0), (-1,-2) and (-1, 2)
c
         include 'DELTAfiles/Ctwo.f'

         W(:,:,:,:,2,it) = Res(:,:,:,:) / 216.0d0
c
c     CASE 3 - quarks at ( 3,0), (-1,-2) and (-1, 2)
c
         include 'DELTAfiles/Cthree.f'

         W(:,:,:,:,3,it) = Res(:,:,:,:) / 216.0d0
c
c     CASE 4 - quarks at ( 3,0), (-2,-3) and (-2, 3)
c
         include 'DELTAfiles/Cfour.f'

         W(:,:,:,:,4,it) = Res(:,:,:,:) / 216.0d0
c
c     CASE 5 - quarks at ( 4,0), (-3,-4) and (-3, 4)
c
c         include 'DELTAfiles/Cfive.f'
c
c         W(:,:,:,:,5,it) = Res(:,:,:,:) / 216.0d0
c
c     CASE 6 - quarks at ( 5,0), (-4,-5) and (-4, 5)
c
c         include 'DELTAfiles/Csix.f'
c
c         W(:,:,:,:,6,it) = Res(:,:,:,:) / 216.0d0
c
c     CASE 7 - quarks at ( 7,0), (-4,-6) and (-4, 6)
c
         include 'DELTAfiles/Cseven.f'

         W(:,:,:,:,5,it) = Res(:,:,:,:) / 216.0d0

c         if(nYloop == 9) then
c
c     CASE 8 - quarks at ( 8,0), (-4,-7) and (-4, 7)
c
c            include 'DELTAfiles/Ceight.f'
c
c            W(:,:,:,:,8,it) = Res(:,:,:,:) / 216.0d0
c
c     CASE 9 - quarks at ( 9,0), (-5,-8) and (-5, 8)
c
c            include 'DELTAfiles/Cnine.f'
c
c            W(:,:,:,:,9,it) = Res(:,:,:,:) / 216.0d0
c
c         endif

      end do

      return

      end subroutine delta_loops

c========================================================================================================

      subroutine delta_loops_mirror(xdir,ydir,ur,ui,Wp,Wm)
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
      include 'DELTAfiles/loop_declarations.f'
c
c     ----------------
c     execution begins
c     ----------------
c
c     -----------------------
c     make the diagonal links
c     -----------------------
c
      call box1x1(xdir,ydir,ur,ui,ld1r,ld1i,rd1r,rd1i)
      call box2x1(xdir,ydir,ur,ui,ld2r,ld2i,rd2r,rd2i)
      call box2x3(xdir,ydir,ld1r,ld1i,rd1r,rd1i,ld2r,ld2i,rd2r,rd2i,ld3r,ld3i,rd3r,rd3i)
c
c     ---------------------------------
c     start of Wilson loop calculations
c     ---------------------------------
c
c     it do loop
c
      do it = 1, nL(that)

         call product(ur,ui,that,tlr,tli,it)
c
c     CASE 1 - quarks at ( 1,0), (-1,-1) and (-1, 1)
c     Mirror - quarks at (-1,0), ( 1, 1) and ( 1,-1)
c
         include 'DELTAfiles/Cone.f'

         Wp(:,:,:,:,1,it) = Res(:,:,:,:) / 216.0d0

         call delta_mirror(xdir,ydir,bxlr,bxli,bllr,blli,brlr,brli,tlr,tli,tl1r,tl1i,Res,it,2,1)

         Wm(:,:,:,:,1,it) = Res(:,:,:,:) / 216.0d0
c
c     CASE 2 - quarks at ( 2,0), (-1,-2) and (-1, 2)
c     Mirror - quarks at (-2,0), ( 1, 2) and ( 1,-2)
c
         include 'DELTAfiles/Ctwo.f'

         Wp(:,:,:,:,2,it) = Res(:,:,:,:) / 216.0d0

         call delta_mirror(xdir,ydir,bxlr,bxli,bllr,blli,brlr,brli,tlr,tli,tl1r,tl1i,Res,it,3,2)

         Wm(:,:,:,:,2,it) = Res(:,:,:,:) / 216.0d0
c
c     CASE 3 - quarks at ( 3,0), (-1,-2) and (-1, 2)
c     Mirror - quarks at (-3,0), ( 1, 2) and ( 1,-2)
c
         include 'DELTAfiles/Cthree.f'

         Wp(:,:,:,:,3,it) = Res(:,:,:,:) / 216.0d0
c
         call delta_mirror(xdir,ydir,bxlr,bxli,bllr,blli,brlr,brli,tlr,tli,tl1r,tl1i,Res,it,4,2)

         Wm(:,:,:,:,3,it) = Res(:,:,:,:) / 216.0d0
c
c     CASE 4 - quarks at ( 3,0), (-2,-3) and (-2, 3)
c     Mirror - quarks at (-3,0), ( 2, 3) and ( 3,-3)
c
         include 'DELTAfiles/Cfour.f'

         Wp(:,:,:,:,4,it) = Res(:,:,:,:) / 216.0d0
c
         call delta_mirror(xdir,ydir,bxlr,bxli,bllr,blli,brlr,brli,tlr,tli,tl1r,tl1i,Res,it,5,3)

         Wm(:,:,:,:,4,it) = Res(:,:,:,:) / 216.0d0
c
c     CASE 5 - quarks at ( 4,0), (-3,-4) and (-3, 4)
c     Mirror - quarks at (-4,0), ( 3, 4) and ( 3,-4)
c
c         include 'DELTAfiles/Cfive.f'
c
c         Wp(:,:,:,:,5,it) = Res(:,:,:,:) / 216.0d0
c
c         call delta_mirror(xdir,ydir,bxlr,bxli,bllr,blli,brlr,brli,tlr,tli,Res,it,4,-3,4)
c
c         Wm(:,:,:,:,5,it) = Res(:,:,:,:) / 216.0d0
c
c     CASE 6 - quarks at ( 5,0), (-4,-5) and (-4, 5)
c     Mirror - quarks at (-5,0), ( 4, 5) and ( 4,-5)
c
c         include 'DELTAfiles/Csix.f'
c
c         Wp(:,:,:,:,6,it) = Res(:,:,:,:) / 216.0d0
c
c         call delta_mirror(xdir,ydir,bxlr,bxli,bllr,blli,brlr,brli,tlr,tli,Res,it,5,-4,5)
c
c         Wm(:,:,:,:,6,it) = Res(:,:,:,:) / 216.0d0
c
c     CASE 7 - quarks at ( 7,0), (-4,-6) and (-4, 6)
c     Mirror - quarks at (-7,0), ( 4, 6) and ( 4,-6)

         include 'DELTAfiles/Cseven.f'

         Wp(:,:,:,:,5,it) = Res(:,:,:,:) / 216.0d0

         call delta_mirror(xdir,ydir,bxlr,bxli,bllr,blli,brlr,brli,tlr,tli,tl1r,tl1i,Res,it,11,6)

         Wm(:,:,:,:,5,it) = Res(:,:,:,:) / 216.0d0
c
c         if(nYloop == 9) then
c
c     CASE 8 - quarks at ( 8,0), (-4,-7) and (-4, 7)
c     Mirror - quarks at (-8,0), ( 4, 7) and ( 4,-7)
c
c            include 'DELTAfiles/Ceight.f'
c
c            Wp(:,:,:,:,8,it) = Res(:,:,:,:) / 216.0d0
c
c            call delta_mirror(xdir,ydir,bxlr,bxli,bllr,blli,brlr,brli,tlr,tli,Res,it,8,-4,7)
c
c            Wm(:,:,:,:,8,it) = Res(:,:,:,:) / 216.0d0
c
c
c     CASE 9 - quarks at ( 9,0), (-5,-8) and (-5, 8)
c     Mirror - quarks at (-9,0), ( 5, 8) and ( 5,-8)
c
c            include 'DELTAfiles/Cnine.f'
c
c            Wp(:,:,:,:,9,it) = Res(:,:,:,:) / 216.0d0
c
c            call delta_mirror(xdir,ydir,bxlr,bxli,bllr,blli,brlr,brli,tlr,tli,Res,it,9,-5,8)
c
c            Wm(:,:,:,:,9,it) = Res(:,:,:,:) / 216.0d0
c
c         endif

      end do

      return

      end subroutine delta_loops_mirror

      END MODULE L_DELTALOOPS
