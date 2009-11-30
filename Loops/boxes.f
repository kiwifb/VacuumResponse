!
!     subroutine to calculate Boxes to represent diagonal links
!     1.0 FB071128
!     2.0 FB090202 added multboxes
!     Derived from the initial code of Adrian Kitson.
!
      MODULE L_BOXES
!
!  Nomeclature:we have 4 kind of boxes mamb mapb pamb papb
!  m is for minus and p is for plus - a and b are integer
!  a is the number of step in x-axis direction
!  b is the number of step in the y direction
!
!----------------------------------------------------
!  We have the four quadrant
!            Top
!        mp   |   pp
! Back       -+-         Front
!        mm   |   pm
!          Bottom
!
!----------------------------------------------------
!
      Contains
!
!  1x1
!
      subroutine boxm1p1(xdir,ydir,ur,ui,boxr,boxi)

      USE L_baryonParam

      IMPLICIT NONE

      include'Templates/Box_dec.f'

      ux1r(:,:,:,:,:,:) = cshift( ur(:,:,:,:,xdir,:,:), dim = xdir, shift =-1)
      ux1i(:,:,:,:,:,:) = cshift( ui(:,:,:,:,xdir,:,:), dim = xdir, shift =-1)

      uy1r(:,:,:,:,:,:) = ur(:,:,:,:,ydir,:,:)
      uy1i(:,:,:,:,:,:) = ui(:,:,:,:,ydir,:,:)

      ux2r(:,:,:,:,:,:) = cshift( ux1r(:,:,:,:,:,:), dim = ydir, shift = 1)
      ux2i(:,:,:,:,:,:) = cshift( ux1i(:,:,:,:,:,:), dim = ydir, shift = 1)

      uy2r(:,:,:,:,:,:) = cshift( uy1r(:,:,:,:,:,:), dim = xdir, shift =-1)
      uy2i(:,:,:,:,:,:) = cshift( uy1i(:,:,:,:,:,:), dim = xdir, shift =-1)

      include'Templates/Box_mp.f'

      return

      end subroutine boxm1p1
!
!----------------------------------------------------
!
      subroutine boxm1m1(xdir,ydir,ur,ui,boxr,boxi)

      USE L_baryonParam

      IMPLICIT NONE

      include'Templates/Box_dec.f'

      ux1r(:,:,:,:,:,:) = cshift( ur(:,:,:,:,xdir,:,:), dim = xdir, shift =-1)
      ux1i(:,:,:,:,:,:) = cshift( ui(:,:,:,:,xdir,:,:), dim = xdir, shift =-1)

      uy1r(:,:,:,:,:,:) = cshift( ur(:,:,:,:,ydir,:,:), dim = ydir, shift =-1)
      uy1i(:,:,:,:,:,:) = cshift( ui(:,:,:,:,ydir,:,:), dim = ydir, shift =-1)

      ux2r(:,:,:,:,:,:) = cshift( ux1r(:,:,:,:,:,:), dim = ydir, shift =-1)
      ux2i(:,:,:,:,:,:) = cshift( ux1i(:,:,:,:,:,:), dim = ydir, shift =-1)

      uy2r(:,:,:,:,:,:) = cshift( uy1r(:,:,:,:,:,:), dim = xdir, shift =-1)
      uy2i(:,:,:,:,:,:) = cshift( uy1i(:,:,:,:,:,:), dim = xdir, shift =-1)

      include'Templates/Box_mm.f'

      return

      end subroutine boxm1m1
!
!----------------------------------------------------
!  2x1
!----------------------------------------------------
!
      subroutine boxm2p1(xdir,ydir,ur,ui,boxr,boxi)

      USE L_baryonParam
      USE L_product

      IMPLICIT NONE

      include'Templates/Box_dec.f'

      call product(ur,ui,xdir,ux1r,ux1i,1)
      call product(ur,ui,xdir,ux1r,ux1i,2)

      ux1r(:,:,:,:,:,:) = cshift( ux1r(:,:,:,:,:,:), dim = xdir, shift =-2)
      ux1i(:,:,:,:,:,:) = cshift( ux1i(:,:,:,:,:,:), dim = xdir, shift =-2)

      uy1r(:,:,:,:,:,:) = ur(:,:,:,:,ydir,:,:)
      uy1i(:,:,:,:,:,:) = ui(:,:,:,:,ydir,:,:)

      ux2r(:,:,:,:,:,:) = cshift( ux1r(:,:,:,:,:,:), dim = ydir, shift = 1)
      ux2i(:,:,:,:,:,:) = cshift( ux1i(:,:,:,:,:,:), dim = ydir, shift = 1)

      uy2r(:,:,:,:,:,:) = cshift( uy1r(:,:,:,:,:,:), dim = xdir, shift =-2)
      uy2i(:,:,:,:,:,:) = cshift( uy1i(:,:,:,:,:,:), dim = xdir, shift =-2)


      include'Templates/Box_mp.f'

      end subroutine boxm2p1
!
!----------------------------------------------------
!
      subroutine boxm2m1(xdir,ydir,ur,ui,boxr,boxi)

      USE L_baryonParam
      USE L_product

      IMPLICIT NONE

      include'Templates/Box_dec.f'

      call product(ur,ui,xdir,ux1r,ux1i,1)
      call product(ur,ui,xdir,ux1r,ux1i,2)

      ux1r(:,:,:,:,:,:) = cshift( ux1r(:,:,:,:,:,:), dim = xdir, shift =-2)
      ux1i(:,:,:,:,:,:) = cshift( ux1i(:,:,:,:,:,:), dim = xdir, shift =-2)

      uy1r(:,:,:,:,:,:) = cshift( ur(:,:,:,:,ydir,:,:), dim = ydir, shift =-1)
      uy1i(:,:,:,:,:,:) = cshift( ui(:,:,:,:,ydir,:,:), dim = ydir, shift =-1)

      ux2r(:,:,:,:,:,:) = cshift( ux1r(:,:,:,:,:,:), dim = ydir, shift =-1)
      ux2i(:,:,:,:,:,:) = cshift( ux1i(:,:,:,:,:,:), dim = ydir, shift =-1)

      uy2r(:,:,:,:,:,:) = cshift( uy1r(:,:,:,:,:,:), dim = xdir, shift =-2)
      uy2i(:,:,:,:,:,:) = cshift( uy1i(:,:,:,:,:,:), dim = xdir, shift =-2)


      include'Templates/Box_mm.f'

      end subroutine boxm2m1
!
!----------------------------------------------------
!  1x2
!----------------------------------------------------
!
      subroutine boxm1p2(xdir,ydir,ur,ui,boxr,boxi)

      USE L_baryonParam
      USE L_product

      IMPLICIT NONE

      include'Templates/Box_dec.f'

      ux1r(:,:,:,:,:,:) = cshift( ur(:,:,:,:,xdir,:,:), dim = xdir, shift =-1)
      ux1i(:,:,:,:,:,:) = cshift( ui(:,:,:,:,xdir,:,:), dim = xdir, shift =-1)

      call product(ur,ui,ydir,uy1r,uy1i,1)
      call product(ur,ui,ydir,uy1r,uy1i,2)

      ux2r(:,:,:,:,:,:) = cshift( ux1r(:,:,:,:,:,:), dim = ydir, shift = 2)
      ux2i(:,:,:,:,:,:) = cshift( ux1i(:,:,:,:,:,:), dim = ydir, shift = 2)

      uy2r(:,:,:,:,:,:) = cshift( uy1r(:,:,:,:,:,:), dim = xdir, shift =-1)
      uy2i(:,:,:,:,:,:) = cshift( uy1i(:,:,:,:,:,:), dim = xdir, shift =-1)


      include'Templates/Box_mp.f'

      end subroutine boxm1p2
!
!----------------------------------------------------
!
      subroutine boxm1m2(xdir,ydir,ur,ui,boxr,boxi)

      USE L_baryonParam
      USE L_product

      IMPLICIT NONE

      include'Templates/Box_dec.f'

      call product(ur,ui,ydir,uy1r,uy1i,1)
      call product(ur,ui,ydir,uy1r,uy1i,2)

      ux1r(:,:,:,:,:,:) = cshift( ur(:,:,:,:,xdir,:,:), dim = xdir, shift =-1)
      ux1i(:,:,:,:,:,:) = cshift( ui(:,:,:,:,xdir,:,:), dim = xdir, shift =-1)

      uy1r(:,:,:,:,:,:) = cshift( uy1r(:,:,:,:,:,:), dim = ydir, shift =-2)
      uy1i(:,:,:,:,:,:) = cshift( uy1i(:,:,:,:,:,:), dim = ydir, shift =-2)

      ux2r(:,:,:,:,:,:) = cshift( ux1r(:,:,:,:,:,:), dim = ydir, shift =-2)
      ux2i(:,:,:,:,:,:) = cshift( ux1i(:,:,:,:,:,:), dim = ydir, shift =-2)

      uy2r(:,:,:,:,:,:) = cshift( uy1r(:,:,:,:,:,:), dim = xdir, shift =-1)
      uy2i(:,:,:,:,:,:) = cshift( uy1i(:,:,:,:,:,:), dim = xdir, shift =-1)


      include'Templates/Box_mm.f'

      end subroutine boxm1m2
!
!----------------------------------------------------
!  3x2
!----------------------------------------------------
!
      subroutine boxm3p2(xdir,ydir,ur,ui,boxr,boxi)

      USE L_baryonParam
      USE L_product

      IMPLICIT NONE

      include'Templates/Box_dec.f'

      call boxm2p1(xdir,ydir,ur,ui,ux1r,ux1i)
      call boxm1p1(xdir,ydir,ur,ui,uy1r,uy1i)

      ux2r(:,:,:,:,:,:) = cshift( cshift( ux1r(:,:,:,:,:,:), dim = ydir, shift = 1), dim = xdir, shift =-1)
      ux2i(:,:,:,:,:,:) = cshift( cshift( ux1i(:,:,:,:,:,:), dim = ydir, shift = 1), dim = xdir, shift =-1)

      uy2r(:,:,:,:,:,:) = cshift( cshift( uy1r(:,:,:,:,:,:), dim = xdir, shift =-2), dim = ydir, shift = 1)
      uy2i(:,:,:,:,:,:) = cshift( cshift( uy1i(:,:,:,:,:,:), dim = xdir, shift =-2), dim = ydir, shift = 1)
!
!     The boxes are already in the right directions, so there is no need for conjugation!
!
      include'Templates/Box_pp.f'

      end subroutine boxm3p2
!
!----------------------------------------------------
!
      subroutine boxm3m2(xdir,ydir,ur,ui,boxr,boxi)

      USE L_baryonParam
      USE L_product

      IMPLICIT NONE

      include'Templates/Box_dec.f'

      call boxm2m1(xdir,ydir,ur,ui,ux1r,ux1i)
      call boxm1m1(xdir,ydir,ur,ui,uy1r,uy1i)

      ux2r(:,:,:,:,:,:) = cshift( cshift( ux1r(:,:,:,:,:,:), dim = ydir, shift =-1), dim = xdir, shift =-1)
      ux2i(:,:,:,:,:,:) = cshift( cshift( ux1i(:,:,:,:,:,:), dim = ydir, shift =-1), dim = xdir, shift =-1)

      uy2r(:,:,:,:,:,:) = cshift( cshift( uy1r(:,:,:,:,:,:), dim = xdir, shift =-2), dim = ydir, shift =-1)
      uy2i(:,:,:,:,:,:) = cshift( cshift( uy1i(:,:,:,:,:,:), dim = xdir, shift =-2), dim = ydir, shift =-1)
!
!     The boxes are already in the right directions, so there is no need for conjugation!
!
      include'Templates/Box_pp.f'

      end subroutine boxm3m2
!
!----------------------------------------------------
!  2x3
!----------------------------------------------------
!
      subroutine boxm2p3(xdir,ydir,ur,ui,boxr,boxi)

      USE L_baryonParam
      USE L_product

      IMPLICIT NONE

      include'Templates/Box_dec.f'

      call boxm1p1(xdir,ydir,ur,ui,ux1r,ux1i)
      call boxm1p2(xdir,ydir,ur,ui,uy1r,uy1i)

      ux2r(:,:,:,:,:,:) = cshift( cshift( ux1r(:,:,:,:,:,:), dim = ydir, shift = 2), dim = xdir, shift =-1)
      ux2i(:,:,:,:,:,:) = cshift( cshift( ux1i(:,:,:,:,:,:), dim = ydir, shift = 2), dim = xdir, shift =-1)

      uy2r(:,:,:,:,:,:) = cshift( cshift( uy1r(:,:,:,:,:,:), dim = xdir, shift =-1), dim = ydir, shift = 1)
      uy2i(:,:,:,:,:,:) = cshift( cshift( uy1i(:,:,:,:,:,:), dim = xdir, shift =-1), dim = ydir, shift = 1)
!
!     The boxes are already in the right directions, so there is no need for conjugation!
!
      include'Templates/Box_pp.f'

      end subroutine boxm2p3
!
!----------------------------------------------------
!
      subroutine boxm2m3(xdir,ydir,ur,ui,boxr,boxi)

      USE L_baryonParam
      USE L_product

      IMPLICIT NONE

      include'Templates/Box_dec.f'

      call boxm1m1(xdir,ydir,ur,ui,ux1r,ux1i)
      call boxm1m2(xdir,ydir,ur,ui,uy1r,uy1i)

      ux2r(:,:,:,:,:,:) = cshift( cshift( ux1r(:,:,:,:,:,:), dim = ydir, shift =-2), dim = xdir, shift =-1)
      ux2i(:,:,:,:,:,:) = cshift( cshift( ux1i(:,:,:,:,:,:), dim = ydir, shift =-2), dim = xdir, shift =-1)

      uy2r(:,:,:,:,:,:) = cshift( cshift( uy1r(:,:,:,:,:,:), dim = xdir, shift =-1), dim = ydir, shift =-1)
      uy2i(:,:,:,:,:,:) = cshift( cshift( uy1i(:,:,:,:,:,:), dim = xdir, shift =-1), dim = ydir, shift =-1)
!
!     The boxes are already in the right directions, so there is no need for conjugation!
!
      include'Templates/Box_pp.f'

      end subroutine boxm2m3
!
! Front links
!
      subroutine boxp1p1(xdir,ydir,ur,ui,boxr,boxi)

      USE L_baryonParam

      IMPLICIT NONE

      include'Templates/Box_dec.f'

      ux1r(:,:,:,:,:,:) = ur(:,:,:,:,xdir,:,:)
      ux1i(:,:,:,:,:,:) = ui(:,:,:,:,xdir,:,:)

      uy1r(:,:,:,:,:,:) = ur(:,:,:,:,ydir,:,:)
      uy1i(:,:,:,:,:,:) = ui(:,:,:,:,ydir,:,:)

      ux2r(:,:,:,:,:,:) = cshift( ux1r(:,:,:,:,:,:), dim = ydir, shift = 1)
      ux2i(:,:,:,:,:,:) = cshift( ux1i(:,:,:,:,:,:), dim = ydir, shift = 1)

      uy2r(:,:,:,:,:,:) = cshift( uy1r(:,:,:,:,:,:), dim = xdir, shift = 1)
      uy2i(:,:,:,:,:,:) = cshift( uy1i(:,:,:,:,:,:), dim = xdir, shift = 1)

      include'Templates/Box_pp.f'

      return

      end subroutine boxp1p1
!
!----------------------------------------------------
!
      subroutine boxp1m1(xdir,ydir,ur,ui,boxr,boxi)

      USE L_baryonParam

      IMPLICIT NONE

      include'Templates/Box_dec.f'

      ux1r(:,:,:,:,:,:) = ur(:,:,:,:,xdir,:,:)
      ux1i(:,:,:,:,:,:) = ui(:,:,:,:,xdir,:,:)

      uy1r(:,:,:,:,:,:) = cshift( ur(:,:,:,:,ydir,:,:), dim = ydir, shift =-1)
      uy1i(:,:,:,:,:,:) = cshift( ui(:,:,:,:,ydir,:,:), dim = ydir, shift =-1)

      ux2r(:,:,:,:,:,:) = cshift( ux1r(:,:,:,:,:,:), dim = ydir, shift =-1)
      ux2i(:,:,:,:,:,:) = cshift( ux1i(:,:,:,:,:,:), dim = ydir, shift =-1)

      uy2r(:,:,:,:,:,:) = cshift( uy1r(:,:,:,:,:,:), dim = xdir, shift = 1)
      uy2i(:,:,:,:,:,:) = cshift( uy1i(:,:,:,:,:,:), dim = xdir, shift = 1)

      include'Templates/Box_pm.f'

      return

      end subroutine boxp1m1
!
!----------------------------------------------------
!
      subroutine boxp2p1(xdir,ydir,ur,ui,boxr,boxi)

      USE L_baryonParam
      USE L_product

      IMPLICIT NONE

      include'Templates/Box_dec.f'

      call product(ur,ui,xdir,ux1r,ux1i,1)
      call product(ur,ui,xdir,ux1r,ux1i,2)

      uy1r(:,:,:,:,:,:) = ur(:,:,:,:,ydir,:,:)
      uy1i(:,:,:,:,:,:) = ui(:,:,:,:,ydir,:,:)

      ux2r(:,:,:,:,:,:) = cshift( ux1r(:,:,:,:,:,:), dim = ydir, shift = 1)
      ux2i(:,:,:,:,:,:) = cshift( ux1i(:,:,:,:,:,:), dim = ydir, shift = 1)

      uy2r(:,:,:,:,:,:) = cshift( uy1r(:,:,:,:,:,:), dim = xdir, shift = 2)
      uy2i(:,:,:,:,:,:) = cshift( uy1i(:,:,:,:,:,:), dim = xdir, shift = 2)


      include'Templates/Box_pp.f'

      end subroutine boxp2p1
!
!----------------------------------------------------
!
      subroutine boxp2m1(xdir,ydir,ur,ui,boxr,boxi)

      USE L_baryonParam
      USE L_product

      IMPLICIT NONE

      include'Templates/Box_dec.f'

      call product(ur,ui,xdir,ux1r,ux1i,1)
      call product(ur,ui,xdir,ux1r,ux1i,2)

      uy1r(:,:,:,:,:,:) = cshift( ur(:,:,:,:,ydir,:,:), dim = ydir, shift =-1)
      uy1i(:,:,:,:,:,:) = cshift( ui(:,:,:,:,ydir,:,:), dim = ydir, shift =-1)

      ux2r(:,:,:,:,:,:) = cshift( ux1r(:,:,:,:,:,:), dim = ydir, shift =-1)
      ux2i(:,:,:,:,:,:) = cshift( ux1i(:,:,:,:,:,:), dim = ydir, shift =-1)

      uy2r(:,:,:,:,:,:) = cshift( uy1r(:,:,:,:,:,:), dim = xdir, shift = 2)
      uy2i(:,:,:,:,:,:) = cshift( uy1i(:,:,:,:,:,:), dim = xdir, shift = 2)


      include'Templates/Box_pm.f'

      end subroutine boxp2m1
!
!----------------------------------------------------
!   General boxes multiplication
!
!  Put box2 at the end of box1 return the result in box1
!  xpos and ypos indicate the end of box1
!
!----------------------------------------------------
!
      subroutine multboxes(xdir,ydir,xpos,ypos,box1r,box1i,box2r,box2i)

      USE L_baryonParam

      IMPLICIT NONE

      integer                                                 :: xpos,ypos
      integer                                                 :: xdir,ydir
      double precision,dimension(nx,ny,nz,nt,nc,nc)           :: box1r,box1i,box2r,box2i
!HPF$ DISTRIBUTE box1r(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE box1i(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE box2r(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE box2i(*,*,BLOCK,BLOCK,*,*)

      integer                                                 ::ic,jc,kc

      double precision,dimension(nx,ny,nz,nt,nc,nc)           :: sbox2r,sbox2i
!HPF$ DISTRIBUTE sbox2r(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE sbox2i(*,*,BLOCK,BLOCK,*,*)
      double precision,dimension(nx,ny,nz,nt,nc,nc)           :: tempr,tempi
!HPF$ DISTRIBUTE tempr(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE tempi(*,*,BLOCK,BLOCK,*,*)

      sbox2r = cshift( cshift( box2r, dim=ydir, shift= ypos), dim=xdir, shift= xpos)
      sbox2i = cshift( cshift( box2i, dim=ydir, shift= ypos), dim=xdir, shift= xpos)

      tempr = 0.d0
      tempi = 0.d0

      do ic = 1,nc
         do jc = 1,nc
            do kc = 1,nc

               tempr(:,:,:,:,ic,jc) = tempr(:,:,:,:,ic,jc)
     &               + box1r(:,:,:,:,ic,kc) * sbox2r(:,:,:,:,kc,jc)
     &               - box1i(:,:,:,:,ic,kc) * sbox2i(:,:,:,:,kc,jc)
               tempi(:,:,:,:,ic,jc) = tempi(:,:,:,:,ic,jc)
     &               + box1r(:,:,:,:,ic,kc) * sbox2i(:,:,:,:,kc,jc)
     &               + box1i(:,:,:,:,ic,kc) * sbox2r(:,:,:,:,kc,jc)

            end do
         end do
      end do

      box1r = tempr
      box1i = tempi

      end subroutine multboxes

      END MODULE L_BOXES
