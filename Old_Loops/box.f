      MODULE L_box
c
c     Create boxes of links for the Y-shape loops
c     Here is the space orientation:
c
c     y ^
c       |
c       o--> x
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
c
      USE L_baryonparam

      Contains

      subroutine box1x1(ur,ui,boxpr,boxpi,boxmr,boxmi)

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

      ux1r(:,:,:,:,:,:) = ur(:,:,:,:,xhat,:,:)
      ux1i(:,:,:,:,:,:) = ui(:,:,:,:,xhat,:,:)

      uy1r(:,:,:,:,:,:) = cshift( ur(:,:,:,:,yhat,:,:), dim = yhat, shift = -1)
      uy1i(:,:,:,:,:,:) = cshift( ui(:,:,:,:,yhat,:,:), dim = yhat, shift = -1)

      ux2r(:,:,:,:,:,:) = cshift( ux1r(:,:,:,:,:,:), dim = yhat, shift = -1)
      ux2i(:,:,:,:,:,:) = cshift( ux1i(:,:,:,:,:,:), dim = yhat, shift = -1)

      uy2r(:,:,:,:,:,:) = cshift( uy1r(:,:,:,:,:,:), dim = xhat, shift = 1)
      uy2i(:,:,:,:,:,:) = cshift( uy1i(:,:,:,:,:,:), dim = xhat, shift = 1)

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
      
      ux1r(:,:,:,:,:,:) = cshift( ux1r(:,:,:,:,:,:), dim = xhat, shift = -1)
      ux1i(:,:,:,:,:,:) = cshift( ux1i(:,:,:,:,:,:), dim = xhat, shift = -1)

      ux2r(:,:,:,:,:,:) = cshift( ux2r(:,:,:,:,:,:), dim = xhat, shift = -1)
      ux2i(:,:,:,:,:,:) = cshift( ux2i(:,:,:,:,:,:), dim = xhat, shift = -1)

      uy2r(:,:,:,:,:,:) = cshift( uy1r(:,:,:,:,:,:), dim = xhat, shift = -1)
      uy2i(:,:,:,:,:,:) = cshift( uy1i(:,:,:,:,:,:), dim = xhat, shift = -1)

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
      ux1r  = 0.d0
      ux1i  = 0.d0

      ux2r(:,:,:,:,:,:) = cshift( ur(:,:,:,:,xhat,:,:), dim = xhat, shift = 1)
      ux2i(:,:,:,:,:,:) = cshift( ui(:,:,:,:,xhat,:,:), dim = xhat, shift = 1)

      do ic = 1, nc
         do jc = 1, nc
            do kc= 1, nc

               ux1r(:,:,:,:,ic,jc) = ux1r(:,:,:,:,ic,jc) + 
     &              ur(:,:,:,:,xhat,ic,kc)*ux2r(:,:,:,:,kc,jc) - ui(:,:,:,:,xhat,ic,kc)*ux2i(:,:,:,:,kc,jc) 

               ux1i(:,:,:,:,ic,jc) = ux1i(:,:,:,:,ic,jc) + 
     &              ur(:,:,:,:,xhat,ic,kc)*ux2i(:,:,:,:,kc,jc) + ui(:,:,:,:,xhat,ic,kc)*ux2r(:,:,:,:,kc,jc) 

            end do
         end do
      end do

      uy1r(:,:,:,:,:,:) = cshift( ur(:,:,:,:,yhat,:,:), dim = yhat, shift = -1)
      uy1i(:,:,:,:,:,:) = cshift( ui(:,:,:,:,yhat,:,:), dim = yhat, shift = -1)

      ux2r(:,:,:,:,:,:) = cshift( ux1r(:,:,:,:,:,:), dim = yhat, shift = -1)
      ux2i(:,:,:,:,:,:) = cshift( ux1i(:,:,:,:,:,:), dim = yhat, shift = -1)

      uy2r(:,:,:,:,:,:) = cshift( uy1r(:,:,:,:,:,:), dim = xhat, shift = 2)
      uy2i(:,:,:,:,:,:) = cshift( uy1i(:,:,:,:,:,:), dim = xhat, shift = 2)

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
      
      ux1r(:,:,:,:,:,:) = cshift( ux1r(:,:,:,:,:,:), dim = xhat, shift = -2)
      ux1i(:,:,:,:,:,:) = cshift( ux1i(:,:,:,:,:,:), dim = xhat, shift = -2)

      ux2r(:,:,:,:,:,:) = cshift( ux2r(:,:,:,:,:,:), dim = xhat, shift = -2)
      ux2i(:,:,:,:,:,:) = cshift( ux2i(:,:,:,:,:,:), dim = xhat, shift = -2)

      uy2r(:,:,:,:,:,:) = cshift( uy1r(:,:,:,:,:,:), dim = xhat, shift = -2)
      uy2i(:,:,:,:,:,:) = cshift( uy1i(:,:,:,:,:,:), dim = xhat, shift = -2)

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
      
      END MODULE L_box
