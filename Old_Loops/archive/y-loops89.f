c
c     subroutine to calculate valules for Wilson loops, W, for 'Y' shape
c     v.3.11 AK040129
c
      MODULE L_YLOOPS89

      integer                                                          :: xdir,ydir
      integer,parameter                                                :: nYloop=9

      CONTAINS

c
c     Create boxes of links for the Y-shape loops
c     Here is the space orientation:
c
c        x ^
c          |
c     y <--o
c
c     box1x1 creates the two following boxes:
c
c     boxp       boxm
c          
c     o->o       o<-o
c     |  |  and  |  |
c     V  V       V  V 
c     o->o       o<-o
c
c     Note that they are not transposed of each other
c
c     box2x1 creates the two following boxes:
c
c      boxp          boxm
c
c     o->o->o       o<-o<-o
c     |     |  and  |     |
c     V     V       V     V 
c     o->o->o       o<-o<-o

c========================================================================================================

      subroutine box1x1(ur,ui,boxpr,boxpi,boxmr,boxmi)

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

      integer                                                 :: ic,jc,kc
      
      boxpr=0.d0
      boxpi=0.d0
      boxmr=0.d0
      boxmi=0.d0

      ux1r(:,:,:,:,:,:) = cshift( ur(:,:,:,:,xdir,:,:), dim = xdir, shift =-1)
      ux1i(:,:,:,:,:,:) = cshift( ui(:,:,:,:,xdir,:,:), dim = xdir, shift =-1)

      uy1r(:,:,:,:,:,:) = ur(:,:,:,:,ydir,:,:)
      uy1i(:,:,:,:,:,:) = ui(:,:,:,:,ydir,:,:)

      ux2r(:,:,:,:,:,:) = cshift( ux1r(:,:,:,:,:,:), dim = ydir, shift = 1)
      ux2i(:,:,:,:,:,:) = cshift( ux1i(:,:,:,:,:,:), dim = ydir, shift = 1)

      uy2r(:,:,:,:,:,:) = cshift( uy1r(:,:,:,:,:,:), dim = xdir, shift =-1)
      uy2i(:,:,:,:,:,:) = cshift( uy1i(:,:,:,:,:,:), dim = xdir, shift =-1)

      do ic = 1, nc
         do jc = 1, nc
            do kc = 1, nc

               boxpr(:,:,:,:,ic,jc) = boxpr(:,:,:,:,ic,jc) + 
     &              ux1r(:,:,:,:,kc,ic)*uy2r(:,:,:,:,kc,jc) + ux1i(:,:,:,:,kc,ic)*uy2i(:,:,:,:,kc,jc) +
     &              uy1r(:,:,:,:,ic,kc)*ux2r(:,:,:,:,jc,kc) + uy1i(:,:,:,:,ic,kc)*ux2i(:,:,:,:,jc,kc)

               boxpi(:,:,:,:,ic,jc) = boxpi(:,:,:,:,ic,jc) + 
     &              ux1r(:,:,:,:,kc,ic)*uy2i(:,:,:,:,kc,jc) - ux1i(:,:,:,:,kc,ic)*uy2r(:,:,:,:,kc,jc) -
     &              uy1r(:,:,:,:,ic,kc)*ux2i(:,:,:,:,jc,kc) + uy1i(:,:,:,:,ic,kc)*ux2r(:,:,:,:,jc,kc)
            end do
         end do
      end do
      
      uy1r(:,:,:,:,:,:) = cshift( ur(:,:,:,:,ydir,:,:), dim = ydir, shift =-1)
      uy1i(:,:,:,:,:,:) = cshift( ui(:,:,:,:,ydir,:,:), dim = ydir, shift =-1)

      ux2r(:,:,:,:,:,:) = cshift( ux1r(:,:,:,:,:,:), dim = ydir, shift =-1)
      ux2i(:,:,:,:,:,:) = cshift( ux1i(:,:,:,:,:,:), dim = ydir, shift =-1)

      uy2r(:,:,:,:,:,:) = cshift( uy1r(:,:,:,:,:,:), dim = xdir, shift =-1)
      uy2i(:,:,:,:,:,:) = cshift( uy1i(:,:,:,:,:,:), dim = xdir, shift =-1)

      do ic = 1, nc
         do jc = 1, nc
            do kc = 1, nc

               boxmr(:,:,:,:,ic,jc) = boxmr(:,:,:,:,ic,jc) + 
     &              ux1r(:,:,:,:,kc,ic)*uy2r(:,:,:,:,jc,kc) - ux1i(:,:,:,:,kc,ic)*uy2i(:,:,:,:,jc,kc) +
     &              uy1r(:,:,:,:,kc,ic)*ux2r(:,:,:,:,jc,kc) - uy1i(:,:,:,:,kc,ic)*ux2i(:,:,:,:,jc,kc)

               boxmi(:,:,:,:,ic,jc) = boxmi(:,:,:,:,ic,jc) - 
     &              ux1r(:,:,:,:,kc,ic)*uy2i(:,:,:,:,jc,kc) - ux1i(:,:,:,:,kc,ic)*uy2r(:,:,:,:,jc,kc) -
     &              uy1r(:,:,:,:,kc,ic)*ux2i(:,:,:,:,jc,kc) - uy1i(:,:,:,:,kc,ic)*ux2r(:,:,:,:,jc,kc)
            end do
         end do
      end do

      return
      
      end subroutine box1x1

c========================================================================================================
      
      subroutine box2x1(ur,ui,boxpr,boxpi,boxmr,boxmi)

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

      integer                                                 :: ic,jc,kc
      
      boxpr = 0.d0
      boxpi = 0.d0
      boxmr = 0.d0
      boxmi = 0.d0
      uy1r  = 0.d0
      uy1i  = 0.d0

      uy2r(:,:,:,:,:,:) = cshift( ur(:,:,:,:,ydir,:,:), dim = ydir, shift = 1)
      uy2i(:,:,:,:,:,:) = cshift( ui(:,:,:,:,ydir,:,:), dim = ydir, shift = 1)

      do ic = 1, nc
         do jc = 1, nc
            do kc= 1, nc

               uy1r(:,:,:,:,ic,jc) = uy1r(:,:,:,:,ic,jc) + 
     &              ur(:,:,:,:,ydir,ic,kc)*uy2r(:,:,:,:,kc,jc) - ui(:,:,:,:,ydir,ic,kc)*uy2i(:,:,:,:,kc,jc) 

               uy1i(:,:,:,:,ic,jc) = uy1i(:,:,:,:,ic,jc) + 
     &              ur(:,:,:,:,ydir,ic,kc)*uy2i(:,:,:,:,kc,jc) + ui(:,:,:,:,ydir,ic,kc)*uy2r(:,:,:,:,kc,jc) 

            end do
         end do
      end do

      ux1r(:,:,:,:,:,:) = cshift( ur(:,:,:,:,xdir,:,:), dim = xdir, shift =-1)
      ux1i(:,:,:,:,:,:) = cshift( ui(:,:,:,:,xdir,:,:), dim = xdir, shift =-1)

      ux2r(:,:,:,:,:,:) = cshift( ux1r(:,:,:,:,:,:), dim = ydir, shift = 2)
      ux2i(:,:,:,:,:,:) = cshift( ux1i(:,:,:,:,:,:), dim = ydir, shift = 2)

      uy2r(:,:,:,:,:,:) = cshift( uy1r(:,:,:,:,:,:), dim = xdir, shift =-1)
      uy2i(:,:,:,:,:,:) = cshift( uy1i(:,:,:,:,:,:), dim = xdir, shift =-1)

      do ic = 1, nc
         do jc = 1, nc
            do kc = 1, nc

               boxpr(:,:,:,:,ic,jc) = boxpr(:,:,:,:,ic,jc) + 
     &              ux1r(:,:,:,:,kc,ic)*uy2r(:,:,:,:,kc,jc) + ux1i(:,:,:,:,kc,ic)*uy2i(:,:,:,:,kc,jc) +
     &              uy1r(:,:,:,:,ic,kc)*ux2r(:,:,:,:,jc,kc) + uy1i(:,:,:,:,ic,kc)*ux2i(:,:,:,:,jc,kc)

               boxpi(:,:,:,:,ic,jc) = boxpi(:,:,:,:,ic,jc) + 
     &              ux1r(:,:,:,:,kc,ic)*uy2i(:,:,:,:,kc,jc) - ux1i(:,:,:,:,kc,ic)*uy2r(:,:,:,:,kc,jc) -
     &              uy1r(:,:,:,:,ic,kc)*ux2i(:,:,:,:,jc,kc) + uy1i(:,:,:,:,ic,kc)*ux2r(:,:,:,:,jc,kc)
            end do
         end do
      end do
      
      uy1r(:,:,:,:,:,:) = cshift( uy1r(:,:,:,:,:,:), dim = ydir, shift =-2)
      uy1i(:,:,:,:,:,:) = cshift( uy1i(:,:,:,:,:,:), dim = ydir, shift =-2)

      ux2r(:,:,:,:,:,:) = cshift( ux1r(:,:,:,:,:,:), dim = ydir, shift =-2)
      ux2i(:,:,:,:,:,:) = cshift( ux1i(:,:,:,:,:,:), dim = ydir, shift =-2)

      uy2r(:,:,:,:,:,:) = cshift( uy1r(:,:,:,:,:,:), dim = xdir, shift =-1)
      uy2i(:,:,:,:,:,:) = cshift( uy1i(:,:,:,:,:,:), dim = xdir, shift =-1)

      do ic = 1, nc
         do jc = 1, nc
            do kc = 1, nc

               boxmr(:,:,:,:,ic,jc) = boxmr(:,:,:,:,ic,jc) + 
     &              ux1r(:,:,:,:,kc,ic)*uy2r(:,:,:,:,jc,kc) - ux1i(:,:,:,:,kc,ic)*uy2i(:,:,:,:,jc,kc) +
     &              uy1r(:,:,:,:,kc,ic)*ux2r(:,:,:,:,jc,kc) - uy1i(:,:,:,:,kc,ic)*ux2i(:,:,:,:,jc,kc)

               boxmi(:,:,:,:,ic,jc) = boxmi(:,:,:,:,ic,jc) - 
     &              ux1r(:,:,:,:,kc,ic)*uy2i(:,:,:,:,jc,kc) - ux1i(:,:,:,:,kc,ic)*uy2r(:,:,:,:,jc,kc) -
     &              uy1r(:,:,:,:,kc,ic)*ux2i(:,:,:,:,jc,kc) - uy1i(:,:,:,:,kc,ic)*ux2r(:,:,:,:,jc,kc)
            end do
         end do
      end do

      return
      
      end subroutine box2x1

c========================================================================================================

      subroutine box2x3(bp1r,bp1i,bm1r,bm1i,bp2r,bp2i,bm2r,bm2i,boxpr,boxpi,boxmr,boxmi)

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

      integer                                                 :: ic,jc,kc
      
      boxpr=0.d0
      boxpi=0.d0

      sb1r = cshift(cshift(bp1r,dim=xdir,shift=-1),dim=ydir,shift= 2)
      sb1i = cshift(cshift(bp1i,dim=xdir,shift=-1),dim=ydir,shift= 2)

      sb2r = cshift(cshift(bp2r,dim=xdir,shift=-1),dim=ydir,shift= 1)
      sb2i = cshift(cshift(bp2i,dim=xdir,shift=-1),dim=ydir,shift= 1)

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

      sb1r = cshift(cshift(bm1r,dim=xdir,shift=-1),dim=ydir,shift=-2)
      sb1i = cshift(cshift(bm1i,dim=xdir,shift=-1),dim=ydir,shift=-2)

      sb2r = cshift(cshift(bm2r,dim=xdir,shift=-1),dim=ydir,shift=-1)
      sb2i = cshift(cshift(bm2i,dim=xdir,shift=-1),dim=ydir,shift=-1)

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

      subroutine y_mirror(bxlr,bxli,bllr,blli,brlr,brli,uptr,upti,Res,it,xtop,xbot,yend)
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

      include 'Yfiles/res.f'            
      
      return

      end subroutine y_mirror
      
c========================================================================================================

      subroutine y_loops_mirror(ur,ui,Wp,Wm)
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
      double precision,dimension(nx,ny,nz,nt,mu,nc,nc)         	 	:: ur,ui
!HPF$ DISTRIBUTE ur(*,*,BLOCK,BLOCK,*,*,*)
!HPF$ DISTRIBUTE ui(*,*,BLOCK,BLOCK,*,*,*)
c
c     the cases (0-nYloop) are stored in the second to last place in 'W'
c
      double precision,dimension(nx,ny,nz,nt,8:nYloop,nL(that))         :: Wp,Wm
!HPF$ DISTRIBUTE Wp(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE Wm(*,*,BLOCK,BLOCK,*,*)
c
      include 'Yfiles/loop_declarations.f'
c
c     ----------------
c     execution begins
c     ----------------
c
c     -----------------------
c     make the diagonal links
c     -----------------------
c
      call box1x1(ur,ui,ld1r,ld1i,rd1r,rd1i)
      call box2x1(ur,ui,ld2r,ld2i,rd2r,rd2i)
      call box2x3(ld1r,ld1i,rd1r,rd1i,ld2r,ld2i,rd2r,rd2i,ld3r,ld3i,rd3r,rd3i)
c
c     ---------------------------------
c     start of Wilson loop calculations
c     ---------------------------------
c
c     it do loop
c
      do it = 1, nL(that)
      
         call product(ur,ui,that,uptr,upti,it)

         call product(ur,ui,xhat,bxlr,bxli,1)
         call product(ur,ui,xhat,bxlr,bxli,2)
         call product(ur,ui,xhat,bxlr,bxli,3)
         call product(ur,ui,xhat,bxlr,bxli,4)
         call product(ur,ui,xhat,bxlr,bxli,5)
         call product(ur,ui,xhat,bxlr,bxli,6)
         call product(ur,ui,xhat,bxlr,bxli,7)
c
c     CASE 8 - quarks at ( 8,0), (-4,-7) and (-4, 7)
c     Mirror - quarks at (-8,0), ( 4, 7) and ( 4,-7)
c     
         include 'Yfiles/Ceight.f'         

         Wp(:,:,:,:,8,it) = Res(:,:,:,:) / 6.0d0
c     
         call y_mirror(bxlr,bxli,bllr,blli,brlr,brli,uptr,upti,Res,it,8,-4,7)

         Wm(:,:,:,:,8,it) = Res(:,:,:,:) / 6.0d0
         
c
c     CASE 9 - quarks at ( 9,0), (-5,-8) and (-5, 8)
c     Mirror - quarks at (-9,0), ( 5, 8) and ( 5,-8)
c     
         include 'Yfiles/Cnine.f'         

         Wp(:,:,:,:,9,it) = Res(:,:,:,:) / 6.0d0
c     
         call y_mirror(bxlr,bxli,bllr,blli,brlr,brli,uptr,upti,Res,it,9,-5,8)

         Wm(:,:,:,:,9,it) = Res(:,:,:,:) / 6.0d0
                     
      end do
      
      return
      
      end subroutine y_loops_mirror

      END MODULE L_YLOOPS89
