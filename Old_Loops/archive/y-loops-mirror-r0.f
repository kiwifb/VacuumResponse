c
c     subroutine to calculate valules for Wilson loops, W, for 'Y' shape
c     v.3.11 AK040129
c
      MODULE YLOOPSMIRROR

      USE baryonParam

      integer                                                          :: xdir,ydir
      integer,parameter                                                :: nYloop=4

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

      subroutine y_loops(ur,ui,Wp,Wm)
c     
c     modules
c     
      USE epsilonIndex
      USE product
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
c     the five cases (0-4) are stored in the second to last place in 'W'
c
      double precision,dimension(nx,ny,nz,nt,0:nYloop,nL(that))         :: Wp,Wm
!HPF$ DISTRIBUTE Wp(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE Wm(*,*,BLOCK,BLOCK,*,*)
c
c     local variables
c
      integer                                                 		:: it
      integer                                                 		:: ic,jc,kc
c
c     diagonal links
c
      double precision,dimension(nx,ny,nz,nt,nc,nc)        		:: rd1r,rd1i,ld1r,ld1i
!HPF$ DISTRIBUTE rd1r(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE rd1i(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE ld1r(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE ld1i(*,*,BLOCK,BLOCK,*,*)
      double precision,dimension(nx,ny,nz,nt,nc,nc)    			:: rd2r,rd2i,ld2r,ld2i
!HPF$ DISTRIBUTE rd2r(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE rd2i(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE ld2r(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE ld2i(*,*,BLOCK,BLOCK,*,*)
c
c     products of links
c
      double precision,dimension(nx,ny,nz,nt,nc,nc)           		:: upxr,upxi
!HPF$ DISTRIBUTE upxr(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE upxi(*,*,BLOCK,BLOCK,*,*)
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
c
c     make products of links 2 units long in the xdir direction
c
      sllr = cshift(ur(:,:,:,:,xdir,:,:),dim=xdir,shift=1)
      slli = cshift(ui(:,:,:,:,xdir,:,:),dim=xdir,shift=1)
      upxr = 0.0d0
      upxi = 0.0d0
      
      do ic = 1,nc
         do jc = 1,nc
            do kc = 1,nc
            
               upxr(:,:,:,:,ic,jc) = upxr(:,:,:,:,ic,jc) 
     &              + ur(:,:,:,:,xdir,ic,kc) * sllr(:,:,:,:,kc,jc) - ui(:,:,:,:,xdir,ic,kc) * slli(:,:,:,:,kc,jc)
             
               upxi(:,:,:,:,ic,jc) = upxi(:,:,:,:,ic,jc)  
     &              + ur(:,:,:,:,xdir,ic,kc) * slli(:,:,:,:,kc,jc) + ui(:,:,:,:,xdir,ic,kc) * sllr(:,:,:,:,kc,jc)
        
            end do
         end do
      end do
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
c     ------------------------------
c     CASE 0 - three quarks at (0,0)
c     ------------------------------
c
         Res = 0.0d0      
         do ic = 1,36

            Res(:,:,:,:) = Res(:,:,:,:)  + fper(ic) * (
     &             uptr(:,:,:,:,la(ic),lap(ic)) * uptr(:,:,:,:,lb(ic),lbp(ic)) * uptr(:,:,:,:,lc(ic),lcp(ic))
     &           - uptr(:,:,:,:,la(ic),lap(ic)) * upti(:,:,:,:,lb(ic),lbp(ic)) * upti(:,:,:,:,lc(ic),lcp(ic))
     &           - upti(:,:,:,:,la(ic),lap(ic)) * uptr(:,:,:,:,lb(ic),lbp(ic)) * upti(:,:,:,:,lc(ic),lcp(ic))
     &           - upti(:,:,:,:,la(ic),lap(ic)) * upti(:,:,:,:,lb(ic),lbp(ic)) * uptr(:,:,:,:,lc(ic),lcp(ic))
     &           )

         end do
         Wp(:,:,:,:,0,it) = Res(:,:,:,:) / 6.0d0
         Wm(:,:,:,:,0,it) = Res(:,:,:,:) / 6.0d0
c
c     ----------------------------------------------
c     CASE 2 - quarks at ( 2,0), (-1,-2) and (-1, 2)
c     Mirror - quarks at (-2,0), ( 1, 2) and ( 1,-2)
c     ----------------------------------------------
c
c     top leg ( 2,0)
c
         tmpr = 0.0d0
         tmpi = 0.0d0
         bxlr(:,:,:,:,:,:) = upxr(:,:,:,:,:,:)
         bxli(:,:,:,:,:,:) = upxi(:,:,:,:,:,:)
         ustr = cshift(uptr(:,:,:,:,:,:),dim=xdir,shift=2) 
         usti = cshift(upti(:,:,:,:,:,:),dim=xdir,shift=2)
         
         include 'Yfiles/topleg.f'
c
c     left leg, i.e., quark at (-1, 2)
c
         tmpr = 0.0d0
         tmpi = 0.0d0
         bllr = ld2r(:,:,:,:,:,:)
         blli = ld2i(:,:,:,:,:,:)
         ustr = cshift(cshift(uptr(:,:,:,:,:,:),dim=ydir,shift= 2),dim=xdir,shift=-1) 
         usti = cshift(cshift(upti(:,:,:,:,:,:),dim=ydir,shift= 2),dim=xdir,shift=-1)
         
         include 'Yfiles/leftleg.f'
c
c     right leg, i.e., quark at (-1,-2)
c
         tmpr = 0.0d0
         tmpi = 0.0d0
         brlr = rd2r(:,:,:,:,:,:)
         brli = rd2i(:,:,:,:,:,:)
         ustr = cshift(cshift(uptr(:,:,:,:,:,:),dim=ydir,shift=-2),dim=xdir,shift=-1) 
         usti = cshift(cshift(upti(:,:,:,:,:,:),dim=ydir,shift=-2),dim=xdir,shift=-1)
         
         include 'Yfiles/rightleg.f'

         include 'Yfiles/res.f'            

         Wp(:,:,:,:,2,it) = Res(:,:,:,:) / 6.0d0
c
c     Mirror top (-2,0)
c
         tmpr = 0.0d0
         tmpi = 0.0d0
         mbxlr(:,:,:,:,:,:) = cshift( bxlr(:,:,:,:,:,:), dim=xdir, shift=-2)
         mbxli(:,:,:,:,:,:) = cshift( bxli(:,:,:,:,:,:), dim=xdir, shift=-2) 

         ustr = cshift(uptr(:,:,:,:,:,:),dim=xdir,shift=-2) 
         usti = cshift(upti(:,:,:,:,:,:),dim=xdir,shift=-2)
         
         include 'Yfiles/Mtopleg.f'
c
c     Mirror left, i.e., quark at ( 1,-2)
c
         tmpr = 0.0d0
         tmpi = 0.0d0
         mbllr(:,:,:,:,:,:) = cshift( cshift( bllr(:,:,:,:,:,:), dim=xdir, shift= 1), dim=ydir, shift=-2)
         mblli(:,:,:,:,:,:) = cshift( cshift( blli(:,:,:,:,:,:), dim=xdir, shift= 1), dim=ydir, shift=-2)

         ustr = cshift(cshift(uptr(:,:,:,:,:,:),dim=ydir,shift=-2),dim=xdir,shift= 1) 
         usti = cshift(cshift(upti(:,:,:,:,:,:),dim=ydir,shift=-2),dim=xdir,shift= 1)
         
         include 'Yfiles/Mleftleg.f'
c
c     Mirror right, i.e., quark at ( 1, 2)
c
         tmpr = 0.0d0
         tmpi = 0.0d0
         mbrlr(:,:,:,:,:,:) = cshift( cshift( brlr(:,:,:,:,:,:), dim=xdir, shift= 1), dim=ydir, shift= 2)
         mbrli(:,:,:,:,:,:) = cshift( cshift( brli(:,:,:,:,:,:), dim=xdir, shift= 1), dim=ydir, shift= 2)

         ustr = cshift(cshift(uptr(:,:,:,:,:,:),dim=ydir,shift= 2),dim=xdir,shift= 1) 
         usti = cshift(cshift(upti(:,:,:,:,:,:),dim=ydir,shift= 2),dim=xdir,shift= 1)
         
         include 'Yfiles/Mrightleg.f'

         include 'Yfiles/res.f'            

         Wm(:,:,:,:,2,it) = Res(:,:,:,:) / 6.0d0
c
c     ----------------------------------------------
c     CASE 1 - quarks at ( 1,0), (-1,-1) and (-1, 1)
c     Mirror - quarks at (-1,0), ( 1, 1) and ( 1,-1)
c     ----------------------------------------------
c
c     top leg ( 1,0)
c
         tmpr = 0.0d0
         tmpi = 0.0d0
         bxlr(:,:,:,:,:,:) = ur(:,:,:,:,xdir,:,:)
         bxli(:,:,:,:,:,:) = ui(:,:,:,:,xdir,:,:)
         ustr = cshift(uptr(:,:,:,:,:,:),dim=xdir,shift=1) 
         usti = cshift(upti(:,:,:,:,:,:),dim=xdir,shift=1)
         
         include 'Yfiles/topleg.f'
c
c     left leg, i.e., quark at (-1, 1)
c
         tmpr = 0.0d0
         tmpi = 0.0d0
         bllr = ld1r(:,:,:,:,:,:)
         blli = ld1i(:,:,:,:,:,:)
         ustr = cshift(cshift(uptr(:,:,:,:,:,:),dim=ydir,shift= 1),dim=xdir,shift=-1) 
         usti = cshift(cshift(upti(:,:,:,:,:,:),dim=ydir,shift= 1),dim=xdir,shift=-1)
         
         include 'Yfiles/leftleg.f'
c
c     right leg, i.e., quark at (-1,-1)
c
         tmpr = 0.0d0
         tmpi = 0.0d0
         brlr = rd1r(:,:,:,:,:,:)
         brli = rd1i(:,:,:,:,:,:)
         ustr = cshift(cshift(uptr(:,:,:,:,:,:),dim=ydir,shift=-1),dim=xdir,shift=-1) 
         usti = cshift(cshift(upti(:,:,:,:,:,:),dim=ydir,shift=-1),dim=xdir,shift=-1)
         
         include 'Yfiles/rightleg.f'

         include 'Yfiles/res.f'            

         Wp(:,:,:,:,1,it) = Res(:,:,:,:) / 6.0d0
c
c     Mirror top (-1,0)
c
         tmpr = 0.0d0
         tmpi = 0.0d0
         mbxlr(:,:,:,:,:,:) = cshift( bxlr(:,:,:,:,:,:), dim=xdir, shift=-1)
         mbxli(:,:,:,:,:,:) = cshift( bxli(:,:,:,:,:,:), dim=xdir, shift=-1)

         ustr = cshift(uptr(:,:,:,:,:,:),dim=xdir,shift=-1) 
         usti = cshift(upti(:,:,:,:,:,:),dim=xdir,shift=-1)
         
         include 'Yfiles/Mtopleg.f'
c
c     Mirror left, i.e., quark at ( 1,-1)
c
         tmpr = 0.0d0
         tmpi = 0.0d0
         mbllr(:,:,:,:,:,:) = cshift( cshift( bllr(:,:,:,:,:,:), dim=xdir, shift= 1), dim=ydir, shift=-1)
         mblli(:,:,:,:,:,:) = cshift( cshift( blli(:,:,:,:,:,:), dim=xdir, shift= 1), dim=ydir, shift=-1)

         ustr = cshift(cshift(uptr(:,:,:,:,:,:),dim=ydir,shift=-1),dim=xdir,shift= 1) 
         usti = cshift(cshift(upti(:,:,:,:,:,:),dim=ydir,shift=-1),dim=xdir,shift= 1)
         
         include 'Yfiles/Mleftleg.f'
c
c     Mirror right, i.e., quark at ( 1, 1)
c
         tmpr = 0.0d0
         tmpi = 0.0d0
         mbrlr(:,:,:,:,:,:) = cshift( cshift( brlr(:,:,:,:,:,:), dim=xdir, shift= 1), dim=ydir, shift= 1)
         mbrli(:,:,:,:,:,:) = cshift( cshift( brli(:,:,:,:,:,:), dim=xdir, shift= 1), dim=ydir, shift= 1)

         ustr = cshift(cshift(uptr(:,:,:,:,:,:),dim=ydir,shift= 1),dim=xdir,shift= 1) 
         usti = cshift(cshift(upti(:,:,:,:,:,:),dim=ydir,shift= 1),dim=xdir,shift= 1)
         
         include 'Yfiles/Mrightleg.f'

         include 'Yfiles/res.f'            

         Wm(:,:,:,:,1,it) = Res(:,:,:,:) / 6.0d0
c
c     ----------------------------------------------
c     CASE 3 - quarks at ( 3,0), (-2,-2) and (-2, 2)
c     Mirror - quarks at (-3,0), ( 2, 2) and ( 2,-2)
c     ----------------------------------------------
c     
c     top leg ( 3,0)
c     
         tmpr = 0.0d0
         tmpi = 0.0d0
         sxlr = cshift(upxr(:,:,:,:,:,:),dim=xdir,shift=1)
         sxli = cshift(upxi(:,:,:,:,:,:),dim=xdir,shift=1)
         
         do ic = 1,nc
            do jc = 1,nc
               do kc = 1,nc
               
                  tmpr(:,:,:,:,ic,jc) = tmpr(:,:,:,:,ic,jc) 
     &                 + bxlr(:,:,:,:,ic,kc) * sxlr(:,:,:,:,kc,jc) - bxli(:,:,:,:,ic,kc) * sxli(:,:,:,:,kc,jc)
                
                  tmpi(:,:,:,:,ic,jc) = tmpi(:,:,:,:,ic,jc)  
     &                 + bxlr(:,:,:,:,ic,kc) * sxli(:,:,:,:,kc,jc) + bxli(:,:,:,:,ic,kc) * sxlr(:,:,:,:,kc,jc)
             
               end do
            end do
         end do
         
         bxlr = tmpr(:,:,:,:,:,:)
         bxli = tmpi(:,:,:,:,:,:)
         tmpr = 0.0d0
         tmpi = 0.0d0
         ustr = cshift(uptr(:,:,:,:,:,:),dim=xdir,shift=3) 
         usti = cshift(upti(:,:,:,:,:,:),dim=xdir,shift=3)
         
         include 'Yfiles/topleg.f'
c
c     left leg, i.e., quark at (-2, 2)
c
         tmpr = 0.0d0
         tmpi = 0.0d0
         sllr = cshift(cshift(ld1r(:,:,:,:,:,:),dim=ydir,shift= 1),dim=xdir,shift=-1)
         slli = cshift(cshift(ld1i(:,:,:,:,:,:),dim=ydir,shift= 1),dim=xdir,shift=-1)
         
         do ic = 1,nc
            do jc = 1,nc
               do kc = 1,nc
               
                  tmpr(:,:,:,:,ic,jc) = tmpr(:,:,:,:,ic,jc) 
     &                 + bllr(:,:,:,:,ic,kc) * sllr(:,:,:,:,kc,jc) - blli(:,:,:,:,ic,kc) * slli(:,:,:,:,kc,jc)
                
                  tmpi(:,:,:,:,ic,jc) = tmpi(:,:,:,:,ic,jc)  
     &                 + bllr(:,:,:,:,ic,kc) * slli(:,:,:,:,kc,jc) + blli(:,:,:,:,ic,kc) * sllr(:,:,:,:,kc,jc)
            
               end do
            end do
         end do
         
         bllr = tmpr(:,:,:,:,:,:)
         blli = tmpi(:,:,:,:,:,:)
         tmpr = 0.0d0
         tmpi = 0.0d0
         ustr = cshift(cshift(uptr(:,:,:,:,:,:),dim=xdir,shift=-2),dim=ydir,shift= 2) 
         usti = cshift(cshift(upti(:,:,:,:,:,:),dim=xdir,shift=-2),dim=ydir,shift= 2)
         
         include 'Yfiles/leftleg.f'
c
c     right leg, i.e., quark at (-2,-2)
c
         tmpr = 0.0d0
         tmpi = 0.0d0
         srlr = cshift(cshift(rd1r(:,:,:,:,:,:),dim=xdir,shift=-1),dim=ydir,shift=-1)
         srli = cshift(cshift(rd1i(:,:,:,:,:,:),dim=xdir,shift=-1),dim=ydir,shift=-1)
         
         do ic = 1,nc
            do jc = 1,nc
               do kc = 1,nc
               
                  tmpr(:,:,:,:,ic,jc) = tmpr(:,:,:,:,ic,jc) 
     &                 + brlr(:,:,:,:,ic,kc) * srlr(:,:,:,:,kc,jc) - brli(:,:,:,:,ic,kc) * srli(:,:,:,:,kc,jc)
                      
                  tmpi(:,:,:,:,ic,jc) = tmpi(:,:,:,:,ic,jc)  
     &                 + brlr(:,:,:,:,ic,kc) * srli(:,:,:,:,kc,jc) + brli(:,:,:,:,ic,kc) * srlr(:,:,:,:,kc,jc)
       
               end do
            end do
         end do
         
         brlr = tmpr(:,:,:,:,:,:)
         brli = tmpi(:,:,:,:,:,:)
         tmpr = 0.0d0
         tmpi = 0.0d0
         ustr = cshift(cshift(uptr(:,:,:,:,:,:),dim=xdir,shift=-2),dim=ydir,shift=-2) 
         usti = cshift(cshift(upti(:,:,:,:,:,:),dim=xdir,shift=-2),dim=ydir,shift=-2)
         
         include 'Yfiles/rightleg.f'
            
         include 'Yfiles/res.f'

         Wp(:,:,:,:,3,it) = Res(:,:,:,:) / 6.0d0
c     
c     Mirror top (-3,0)
c     
         tmpr = 0.0d0
         tmpi = 0.0d0
         mbxlr(:,:,:,:,:,:) = cshift( bxlr(:,:,:,:,:,:), dim=xdir, shift=-3)
         mbxli(:,:,:,:,:,:) = cshift( bxli(:,:,:,:,:,:), dim=xdir, shift=-3)

         ustr = cshift(uptr(:,:,:,:,:,:),dim=xdir,shift=-3) 
         usti = cshift(upti(:,:,:,:,:,:),dim=xdir,shift=-3)
         
         include 'Yfiles/Mtopleg.f'
c
c     Mirror left, i.e., quark at ( 2,-2)
c
         tmpr = 0.0d0
         tmpi = 0.0d0
         mbllr(:,:,:,:,:,:) = cshift( cshift( bllr(:,:,:,:,:,:), dim=xdir, shift= 2), dim=ydir, shift=-2)
         mblli(:,:,:,:,:,:) = cshift( cshift( blli(:,:,:,:,:,:), dim=xdir, shift= 2), dim=ydir, shift=-2)

         ustr = cshift(cshift(uptr(:,:,:,:,:,:),dim=xdir,shift= 2),dim=ydir,shift=-2) 
         usti = cshift(cshift(upti(:,:,:,:,:,:),dim=xdir,shift= 2),dim=ydir,shift=-2)
         
         include 'Yfiles/Mleftleg.f'
c
c     Mirror right, i.e., quark at ( 2, 2)
c
         tmpr = 0.0d0
         tmpi = 0.0d0
         mbrlr(:,:,:,:,:,:) = cshift( cshift( brlr(:,:,:,:,:,:), dim=xdir, shift= 2), dim=ydir, shift= 2)
         mblli(:,:,:,:,:,:) = cshift( cshift( brli(:,:,:,:,:,:), dim=xdir, shift= 2), dim=ydir, shift= 2)

         ustr = cshift(cshift(uptr(:,:,:,:,:,:),dim=xdir,shift= 2),dim=ydir,shift= 2) 
         usti = cshift(cshift(upti(:,:,:,:,:,:),dim=xdir,shift= 2),dim=ydir,shift= 2)
         
         include 'Yfiles/Mrightleg.f'
            
         include 'Yfiles/res.f'

         Wm(:,:,:,:,3,it) = Res(:,:,:,:) / 6.0d0
c
c     ----------------------------------------------
c     CASE 4 - quarks at ( 5,0), (-3,-4) and (-3, 4)
c     Mirror - quarks at (-5,0), ( 3, 4) and ( 3,-4)
c     ----------------------------------------------
c     
c     top leg ( 5,0)
c     
         tmpr = 0.0d0
         tmpi = 0.0d0
         sxlr = cshift(upxr(:,:,:,:,:,:),dim=xdir,shift=3)
         sxli = cshift(upxi(:,:,:,:,:,:),dim=xdir,shift=3)
         
         do ic = 1,nc
            do jc = 1,nc
               do kc = 1,nc
               
                  tmpr(:,:,:,:,ic,jc) = tmpr(:,:,:,:,ic,jc) 
     &                 + bxlr(:,:,:,:,ic,kc) * sxlr(:,:,:,:,kc,jc) - bxli(:,:,:,:,ic,kc) * sxli(:,:,:,:,kc,jc)
           
                  tmpi(:,:,:,:,ic,jc) = tmpi(:,:,:,:,ic,jc)  
     &                 + bxlr(:,:,:,:,ic,kc) * sxli(:,:,:,:,kc,jc) + bxli(:,:,:,:,ic,kc) * sxlr(:,:,:,:,kc,jc)
             
               end do
            end do
         end do
         
         bxlr = tmpr(:,:,:,:,:,:)
         bxli = tmpi(:,:,:,:,:,:)
         tmpr = 0.0d0
         tmpi = 0.0d0
         ustr = cshift(uptr(:,:,:,:,:,:),dim=ydir,shift=5) 
         usti = cshift(upti(:,:,:,:,:,:),dim=ydir,shift=5)
         
         include 'Yfiles/topleg.f'
c
c     left leg, i.e., quark at (-3, 4)
c
         tmpr = 0.0d0
         tmpi = 0.0d0
         sllr = cshift(cshift(ld2r(:,:,:,:,:,:),dim=xdir,shift=-2),dim=ydir,shift= 2)
         slli = cshift(cshift(ld2i(:,:,:,:,:,:),dim=xdir,shift=-2),dim=ydir,shift= 2)
         
         do ic = 1,nc
            do jc = 1,nc
               do kc = 1,nc
                  
                  tmpr(:,:,:,:,ic,jc) = tmpr(:,:,:,:,ic,jc) 
     &                 + bllr(:,:,:,:,ic,kc) * sllr(:,:,:,:,kc,jc) - blli(:,:,:,:,ic,kc) * slli(:,:,:,:,kc,jc)
       
                  tmpi(:,:,:,:,ic,jc) = tmpi(:,:,:,:,ic,jc)  
     &                 + bllr(:,:,:,:,ic,kc) * slli(:,:,:,:,kc,jc) + blli(:,:,:,:,ic,kc) * sllr(:,:,:,:,kc,jc)
              
               end do
            end do
         end do
         
         bllr = tmpr(:,:,:,:,:,:)
         blli = tmpi(:,:,:,:,:,:)
         tmpr = 0.0d0
         tmpi = 0.0d0
         ustr = cshift(cshift(uptr(:,:,:,:,:,:),dim=xdir,shift=-3),dim=ydir,shift= 4) 
         usti = cshift(cshift(upti(:,:,:,:,:,:),dim=xdir,shift=-3),dim=ydir,shift= 4)
         
         include 'Yfiles/leftleg.f'
c     
c     right leg, i.e., quark at (-3,-4)
c
         tmpr = 0.0d0
         tmpi = 0.0d0
         srlr = cshift(cshift(rd2r(:,:,:,:,:,:),dim=xdir,shift=-2),dim=ydir,shift=-2)
         srli = cshift(cshift(rd2i(:,:,:,:,:,:),dim=xdir,shift=-2),dim=ydir,shift=-2)
         
         do ic = 1,nc
            do jc = 1,nc
               do kc = 1,nc
               
                  tmpr(:,:,:,:,ic,jc) = tmpr(:,:,:,:,ic,jc) 
     &                 + brlr(:,:,:,:,ic,kc) * srlr(:,:,:,:,kc,jc) - brli(:,:,:,:,ic,kc) * srli(:,:,:,:,kc,jc)
       
                  tmpi(:,:,:,:,ic,jc) = tmpi(:,:,:,:,ic,jc)  
     &                 + brlr(:,:,:,:,ic,kc) * srli(:,:,:,:,kc,jc) + brli(:,:,:,:,ic,kc) * srlr(:,:,:,:,kc,jc)
            
               end do
            end do
         end do
         
         brlr = tmpr(:,:,:,:,:,:)
         brli = tmpi(:,:,:,:,:,:)
         tmpr = 0.0d0
         tmpi = 0.0d0
         ustr = cshift(cshift(uptr(:,:,:,:,:,:),dim=xdir,shift=-3),dim=ydir,shift=-4) 
         usti = cshift(cshift(upti(:,:,:,:,:,:),dim=xdir,shift=-3),dim=ydir,shift=-4)
         
         include 'Yfiles/rightleg.f'

         include 'Yfiles/res.f'         

         Wp(:,:,:,:,4,it) = Res(:,:,:,:) / 6.0d0
c     
c     Mirror top (-5,0)
c     
         tmpr = 0.0d0
         tmpi = 0.0d0
         mbxlr(:,:,:,:,:,:) = cshift( bxlr(:,:,:,:,:,:), dim=xdir, shift=-5)
         mbxli(:,:,:,:,:,:) = cshift( bxli(:,:,:,:,:,:), dim=xdir, shift=-5)

         ustr = cshift(uptr(:,:,:,:,:,:),dim=ydir,shift=-5) 
         usti = cshift(upti(:,:,:,:,:,:),dim=ydir,shift=-5)
         
         include 'Yfiles/Mtopleg.f'
c
c     Mirror left, i.e., quark at ( 3,-4)
c
         tmpr = 0.0d0
         tmpi = 0.0d0
         mbllr(:,:,:,:,:,:) = cshift( cshift( bllr(:,:,:,:,:,:), dim=xdir, shift= 3), dim=ydir, shift=-4)
         mblli(:,:,:,:,:,:) = cshift( cshift( blli(:,:,:,:,:,:), dim=xdir, shift= 3), dim=ydir, shift=-4)

         ustr = cshift(cshift(uptr(:,:,:,:,:,:),dim=xdir,shift= 3),dim=ydir,shift=-4) 
         usti = cshift(cshift(upti(:,:,:,:,:,:),dim=xdir,shift= 3),dim=ydir,shift=-4)
         
         include 'Yfiles/Mleftleg.f'
c     
c     Mirror right, i.e., quark at ( 3, 4)
c
         tmpr = 0.0d0
         tmpi = 0.0d0
         mbrlr(:,:,:,:,:,:) = cshift( cshift( brlr(:,:,:,:,:,:), dim=xdir, shift= 3), dim=ydir, shift= 4)
         mblli(:,:,:,:,:,:) = cshift( cshift( brli(:,:,:,:,:,:), dim=xdir, shift= 3), dim=ydir, shift= 4)

         ustr = cshift(cshift(uptr(:,:,:,:,:,:),dim=xdir,shift= 3),dim=ydir,shift= 4) 
         usti = cshift(cshift(upti(:,:,:,:,:,:),dim=xdir,shift= 3),dim=ydir,shift= 4)
         
         include 'Yfiles/Mrightleg.f'

         include 'Yfiles/res.f'         

         Wm(:,:,:,:,4,it) = Res(:,:,:,:) / 6.0d0
      
      end do
      
      return
      
      end subroutine y_loops

      END MODULE YLOOPSMIRROR
