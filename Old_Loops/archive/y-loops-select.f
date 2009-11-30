c
c     subroutine to calculate valules for Wilson loops, W, for 'Y' shape
c     v.3.11 AK040129
c
      MODULE YLOOPSSELECT

      CONTAINS

      subroutine y_loops(ur,ui,W)
c     
c     modules
c     
      USE baryonParam
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
      double precision,dimension(nx,ny,nz,nt,0:4,1:nL(that)) 		:: W
!HPF$ DISTRIBUTE W(*,*,BLOCK,BLOCK,*,*)
c
c     local variables
c
      integer                                                 		:: it,ib,id
      integer                                                 		:: ic,jc,kc
      integer,dimension(0:255,8)                                        :: BD                                         
c
c     diagonal links
c
      double precision,dimension(nx,ny,nz,nt,0:1,nc,nc)        		:: rd1r,rd1i,ld1r,ld1i
!HPF$ DISTRIBUTE rd1r(*,*,BLOCK,BLOCK,*,*,*)
!HPF$ DISTRIBUTE rd1i(*,*,BLOCK,BLOCK,*,*,*)
!HPF$ DISTRIBUTE ld1r(*,*,BLOCK,BLOCK,*,*,*)
!HPF$ DISTRIBUTE ld1i(*,*,BLOCK,BLOCK,*,*,*)
      double precision,dimension(nx,ny,nz,nt,0:1,nc,nc)    		:: rd2r,rd2i,ld2r,ld2i
!HPF$ DISTRIBUTE rd2r(*,*,BLOCK,BLOCK,*,*,*)
!HPF$ DISTRIBUTE rd2i(*,*,BLOCK,BLOCK,*,*,*)
!HPF$ DISTRIBUTE ld2r(*,*,BLOCK,BLOCK,*,*,*)
!HPF$ DISTRIBUTE ld2i(*,*,BLOCK,BLOCK,*,*,*)
c
c     products of links
c
      double precision,dimension(nx,ny,nz,nt,nc,nc)           		:: upxr,upxi
!HPF$ DISTRIBUTE upxr(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE upxi(*,*,BLOCK,BLOCK,*,*)
      double precision,dimension(nx,ny,nz,nt,nc,nc)           		:: upyr,upyi
!HPF$ DISTRIBUTE upyr(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE upyi(*,*,BLOCK,BLOCK,*,*)
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
c     shifted time links
c
      double precision,dimension(nx,ny,nz,nt,nc,nc)           		:: ustr,usti
!HPF$ DISTRIBUTE ustr(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE usti(*,*,BLOCK,BLOCK,*,*)
c
c     staples
c
      double precision,dimension(nx,ny,nz,nt,nc,nc)        		:: sxlr,sxli
!HPF$ DISTRIBUTE sxlr(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE sxli(*,*,BLOCK,BLOCK,*,*)
      double precision,dimension(nx,ny,nz,nt,nc,nc)        		:: sllr,slli
!HPF$ DISTRIBUTE sllr(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE slli(*,*,BLOCK,BLOCK,*,*)
      double precision,dimension(nx,ny,nz,nt,nc,nc)             	:: srlr,srli
!HPF$ DISTRIBUTE srlr(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE srli(*,*,BLOCK,BLOCK,*,*)
c
c     temporary arrays
c
      double precision,dimension(nx,ny,nz,nt,nc,nc)           		:: tmp1r,tmp1i,tmp2r,tmp2i,tmp3r,tmp3i
!HPF$ DISTRIBUTE tmp1r(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE tmp1i(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE tmp2r(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE tmp2i(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE tmp3r(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE tmp3i(*,*,BLOCK,BLOCK,*,*)
      double precision,dimension(nx,ny,nz,nt)                           :: Res
!HPF$ DISTRIBUTE Res(*,*,BLOCK,BLOCK)
c
c     ----------------
c     execution begins
c     ----------------
c
c     permutations
c
      call setEpsilonIndex()        
c
c     initialize BD
c
      do ic = 0, 255

         ib = ic

         do id = 1, 4

            BD(ic,id)   = mod(ib,2)
            ib = (ib - BD(ic,id))/2

            BD(ic,id+4) = mod(ib,2)
            ib = (ib - BD(ic,id+4))/2
            
         end do

         write(*,'(i3,a,8(i1))') ic,' = ',BD(ic,:)

      end do
c
c     -----------------------
c     make the diagonal links
c     -----------------------
c
      rd1r = 0.0d0
      rd1i = 0.0d0
      ld1r = 0.0d0
      ld1i = 0.0d0
      rd2r = 0.0d0
      rd2i = 0.0d0
      ld2r = 0.0d0
      ld2i = 0.0d0
c     
c     temporary shifted links (cshifts OUTSIDE do loops!)
c
c     the labeling here is confusing as sxl, sll, etc are the 
c     labels for the staples, however here they're just being 
c     used as temporary arrays
c
      sxlr(:,:,:,:,:,:) = cshift(ur(:,:,:,:,yhat,:,:),dim=yhat,shift=-1) 
      sxli(:,:,:,:,:,:) = cshift(ui(:,:,:,:,yhat,:,:),dim=yhat,shift=-1)
      bxlr = cshift(ur(:,:,:,:,xhat,:,:),dim=xhat,shift=-1)
      bxli = cshift(ui(:,:,:,:,xhat,:,:),dim=xhat,shift=-1)
      tmp2r(:,:,:,:,:,:) = cshift(sxlr(:,:,:,:,:,:),dim=xhat,shift=-1)
      tmp2i(:,:,:,:,:,:) = cshift(sxli(:,:,:,:,:,:),dim=xhat,shift=-1)
      tmp3r(:,:,:,:,:,:) = cshift(sxlr(:,:,:,:,:,:),dim=xhat,shift=1)
      tmp3i(:,:,:,:,:,:) = cshift(sxli(:,:,:,:,:,:),dim=xhat,shift=1)
      tmp1r = cshift(ur(:,:,:,:,xhat,:,:),dim=yhat,shift=-1)
      tmp1i = cshift(ui(:,:,:,:,xhat,:,:),dim=yhat,shift=-1)
      bllr(:,:,:,:,:,:) = cshift(bxlr(:,:,:,:,:,:),dim=yhat,shift=-1)
      blli(:,:,:,:,:,:) = cshift(bxli(:,:,:,:,:,:),dim=yhat,shift=-1)
c
c     1x1 square boxes
c
      do ic = 1,nc
         do jc = 1,nc
            do kc = 1,nc
c
c     ld1 left diagonal 1x1
c
               ld1r(:,:,:,:,0,ic,jc) = ld1r(:,:,:,:,0,ic,jc) 
     &              + bxlr(:,:,:,:,kc,ic) * tmp2r(:,:,:,:,jc,kc) 
     &              - bxli(:,:,:,:,kc,ic) * tmp2i(:,:,:,:,jc,kc)
               ld1r(:,:,:,:,1,ic,jc) = ld1r(:,:,:,:,1,ic,jc) 
     &              + sxlr(:,:,:,:,kc,ic) * bllr(:,:,:,:,jc,kc) 
     &              - sxli(:,:,:,:,kc,ic) * blli(:,:,:,:,jc,kc)

               ld1i(:,:,:,:,0,ic,jc) = ld1i(:,:,:,:,0,ic,jc)  
     &              - bxlr(:,:,:,:,kc,ic) * tmp2i(:,:,:,:,jc,kc) 
     &              - bxli(:,:,:,:,kc,ic) * tmp2r(:,:,:,:,jc,kc)
               ld1i(:,:,:,:,1,ic,jc) = ld1i(:,:,:,:,1,ic,jc)  
     &              - sxlr(:,:,:,:,kc,ic) * blli(:,:,:,:,jc,kc) 
     &              - sxli(:,:,:,:,kc,ic) * bllr(:,:,:,:,jc,kc)
c     
c     rd1 right diagonal 1x1
c     
               rd1r(:,:,:,:,0,ic,jc) = rd1r(:,:,:,:,0,ic,jc) 
     &              + ur(:,:,:,:,xhat,ic,kc) * tmp3r(:,:,:,:,jc,kc) 
     &              + ui(:,:,:,:,xhat,ic,kc) * tmp3i(:,:,:,:,jc,kc)
               rd1r(:,:,:,:,1,ic,jc) = rd1r(:,:,:,:,1,ic,jc) 
     &              + sxlr(:,:,:,:,kc,ic) * tmp1r(:,:,:,:,kc,jc) 
     &              + sxli(:,:,:,:,kc,ic) * tmp1i(:,:,:,:,kc,jc)

               rd1i(:,:,:,:,0,ic,jc) = rd1i(:,:,:,:,0,ic,jc)  
     &              - ur(:,:,:,:,xhat,ic,kc) * tmp3i(:,:,:,:,jc,kc) 
     &              + ui(:,:,:,:,xhat,ic,kc) * tmp3r(:,:,:,:,jc,kc)
               rd1i(:,:,:,:,1,ic,jc) = rd1i(:,:,:,:,1,ic,jc)  
     &              + sxlr(:,:,:,:,kc,ic) * tmp1i(:,:,:,:,kc,jc) 
     &              - sxli(:,:,:,:,kc,ic) * tmp1r(:,:,:,:,kc,jc)
             
            end do
         end do
      end do
c
c     first make products of links 2 units long in the xhat direction
c
      tmp2r(:,:,:,:,:,:) = cshift(ur(:,:,:,:,xhat,:,:),dim=xhat,shift=1)
      tmp2i(:,:,:,:,:,:) = cshift(ui(:,:,:,:,xhat,:,:),dim=xhat,shift=1)
      upxr = 0.0d0
      upxi = 0.0d0
      tmp3r(:,:,:,:,:,:) = cshift(ur(:,:,:,:,yhat,:,:),dim=yhat,shift=1)
      tmp3i(:,:,:,:,:,:) = cshift(ui(:,:,:,:,yhat,:,:),dim=yhat,shift=1)
      upyr = 0.0d0
      upyi = 0.0d0
      
      do ic = 1,nc
         do jc = 1,nc
            do kc = 1,nc
            
               upxr(:,:,:,:,ic,jc) = upxr(:,:,:,:,ic,jc) 
     &              + ur(:,:,:,:,xhat,ic,kc) * tmp2r(:,:,:,:,kc,jc) 
     &              - ui(:,:,:,:,xhat,ic,kc) * tmp2i(:,:,:,:,kc,jc)
             
               upxi(:,:,:,:,ic,jc) = upxi(:,:,:,:,ic,jc)  
     &              + ur(:,:,:,:,xhat,ic,kc) * tmp2i(:,:,:,:,kc,jc) 
     &              + ui(:,:,:,:,xhat,ic,kc) * tmp2r(:,:,:,:,kc,jc)
c
c     ...and yhat direction too (used later)
c
               upyr(:,:,:,:,ic,jc) = upyr(:,:,:,:,ic,jc)
     &              + ur(:,:,:,:,yhat,ic,kc) * tmp3r(:,:,:,:,kc,jc) 
     &              - ui(:,:,:,:,yhat,ic,kc) * tmp3i(:,:,:,:,kc,jc)
           
               upyi(:,:,:,:,ic,jc) = upyi(:,:,:,:,ic,jc)  
     &              + ur(:,:,:,:,yhat,ic,kc) * tmp3i(:,:,:,:,kc,jc) 
     &              + ui(:,:,:,:,yhat,ic,kc) * tmp3r(:,:,:,:,kc,jc)
        
            end do
         end do
      end do
c
c     2x1 rectangle boxes
c
      bxlr = cshift(upxr(:,:,:,:,:,:),dim=xhat,shift=-2)
      bxli = cshift(upxi(:,:,:,:,:,:),dim=xhat,shift=-2)
      tmp2r(:,:,:,:,:,:) = cshift(sxlr(:,:,:,:,:,:),dim=xhat,shift=-2)
      tmp2i(:,:,:,:,:,:) = cshift(sxli(:,:,:,:,:,:),dim=xhat,shift=-2)
      tmp3r(:,:,:,:,:,:) = cshift(sxlr(:,:,:,:,:,:),dim=xhat,shift=2)
      tmp3i(:,:,:,:,:,:) = cshift(sxli(:,:,:,:,:,:),dim=xhat,shift=2)
      tmp1r = cshift(upxr(:,:,:,:,:,:),dim=yhat,shift=-1)
      tmp1i = cshift(upxi(:,:,:,:,:,:),dim=yhat,shift=-1)
      bllr(:,:,:,:,:,:) = cshift(bxlr(:,:,:,:,:,:),dim=yhat,shift=-1)
      blli(:,:,:,:,:,:) = cshift(bxli(:,:,:,:,:,:),dim=yhat,shift=-1)
      
      do ic = 1,nc
         do jc = 1,nc
            do kc = 1,nc
c
c     ld2 left diagonal 2x1
c
               ld2r(:,:,:,:,0,ic,jc) = ld2r(:,:,:,:,0,ic,jc) 
     &              + bxlr(:,:,:,:,kc,ic) * tmp2r(:,:,:,:,jc,kc) 
     &              - bxli(:,:,:,:,kc,ic) * tmp2i(:,:,:,:,jc,kc)
               ld2r(:,:,:,:,1,ic,jc) = ld2r(:,:,:,:,1,ic,jc) 
     &              + sxlr(:,:,:,:,kc,ic) * bllr(:,:,:,:,jc,kc) 
     &              - sxli(:,:,:,:,kc,ic) * blli(:,:,:,:,jc,kc)

               ld2i(:,:,:,:,0,ic,jc) = ld2i(:,:,:,:,0,ic,jc)  
     &              - bxlr(:,:,:,:,kc,ic) * tmp2i(:,:,:,:,jc,kc) 
     &              - bxli(:,:,:,:,kc,ic) * tmp2r(:,:,:,:,jc,kc)
               ld2i(:,:,:,:,1,ic,jc) = ld2i(:,:,:,:,1,ic,jc)  
     &              - sxlr(:,:,:,:,kc,ic) * blli(:,:,:,:,jc,kc) 
     &              - sxli(:,:,:,:,kc,ic) * bllr(:,:,:,:,jc,kc)

c
c     rd2 right diagonal 2x1
c
               rd2r(:,:,:,:,0,ic,jc) = rd2r(:,:,:,:,0,ic,jc) 
     &              + upxr(:,:,:,:,ic,kc) * tmp3r(:,:,:,:,jc,kc) 
     &              + upxi(:,:,:,:,ic,kc) * tmp3i(:,:,:,:,jc,kc)
               rd2r(:,:,:,:,1,ic,jc) = rd2r(:,:,:,:,1,ic,jc) 
     &              + sxlr(:,:,:,:,kc,ic) * tmp1r(:,:,:,:,kc,jc) 
     &              + sxli(:,:,:,:,kc,ic) * tmp1i(:,:,:,:,kc,jc)
                   
               rd2i(:,:,:,:,0,ic,jc) = rd2i(:,:,:,:,0,ic,jc)  
     &              - upxr(:,:,:,:,ic,kc) * tmp3i(:,:,:,:,jc,kc) 
     &              + upxi(:,:,:,:,ic,kc) * tmp3r(:,:,:,:,jc,kc)
               rd2i(:,:,:,:,1,ic,jc) = rd2i(:,:,:,:,1,ic,jc)  
     &              + sxlr(:,:,:,:,kc,ic) * tmp1i(:,:,:,:,kc,jc) 
     &              - sxli(:,:,:,:,kc,ic) * tmp1r(:,:,:,:,kc,jc)
  
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
         W(:,:,:,:,0,it) = Res(:,:,:,:) / 6.0d0

         W(:,:,:,:,1:4,it) = 0.0d0

         do ib = 0, 255
c         ib=0

            bllr = ld1r(:,:,:,:,BD(ib,1),:,:)
            blli = ld1i(:,:,:,:,BD(ib,1),:,:)
            brlr = rd1r(:,:,:,:,BD(ib,5),:,:)
            brli = rd1i(:,:,:,:,BD(ib,5),:,:)

            if(ib < 4) then
c
c     --------------------------------------------
c     CASE 1 - quarks at (0,1), (-1,-1) and (1,-1)
c     --------------------------------------------
c
c     top leg
c
               tmp1r = 0.0d0
               tmp1i = 0.0d0
               bxlr(:,:,:,:,:,:) = ur(:,:,:,:,yhat,:,:)
               bxli(:,:,:,:,:,:) = ui(:,:,:,:,yhat,:,:)
               ustr = cshift(uptr(:,:,:,:,:,:),dim=yhat,shift=1) 
               usti = cshift(upti(:,:,:,:,:,:),dim=yhat,shift=1)
         
               include 'Yfiles/topleg2.f'
c
c     left leg, i.e., quark at (-1,-1)
c
               tmp1r = 0.0d0
               tmp1i = 0.0d0
               ustr = cshift(cshift(uptr(:,:,:,:,:,:),dim=xhat,shift=-1),dim=yhat,shift=-1) 
               usti = cshift(cshift(upti(:,:,:,:,:,:),dim=xhat,shift=-1),dim=yhat,shift=-1)
         
               include 'Yfiles/leftleg2.f'
c
c     right leg, i.e., quark at (1,-1)
c
               tmp1r = 0.0d0
               tmp1i = 0.0d0
               ustr = cshift(cshift(uptr(:,:,:,:,:,:),dim=xhat,shift=1),dim=yhat,shift=-1) 
               usti = cshift(cshift(upti(:,:,:,:,:,:),dim=xhat,shift=1),dim=yhat,shift=-1)
         
               include 'Yfiles/rightleg2.f'

               include 'Yfiles/res.f'            

               W(:,:,:,:,1,it) = W(:,:,:,:,1,it) + Res(:,:,:,:) / 6.0d0

            endif

            tmp2r = cshift(cshift(ld1r(:,:,:,:,BD(ib,2),:,:),dim=xhat,shift=-1),dim=yhat,shift=-1)
            tmp2i = cshift(cshift(ld1i(:,:,:,:,BD(ib,2),:,:),dim=xhat,shift=-1),dim=yhat,shift=-1)
            tmp3r = cshift(cshift(rd1r(:,:,:,:,BD(ib,6),:,:),dim=xhat,shift= 1),dim=yhat,shift=-1)
            tmp3i = cshift(cshift(rd1i(:,:,:,:,BD(ib,6),:,:),dim=xhat,shift= 1),dim=yhat,shift=-1)
            tmp1r = 0.0d0
            tmp1i = 0.0d0
            do ic = 1,nc
               do jc = 1,nc
                  do kc = 1,nc
                     
                     tmp1r(:,:,:,:,ic,jc) = tmp1r(:,:,:,:,ic,jc) 
     &                    + bllr(:,:,:,:,ic,kc) * tmp2r(:,:,:,:,kc,jc) 
     &                    - blli(:,:,:,:,ic,kc) * tmp2i(:,:,:,:,kc,jc)
                
                     tmp1i(:,:,:,:,ic,jc) = tmp1i(:,:,:,:,ic,jc)  
     &                    + bllr(:,:,:,:,ic,kc) * tmp2i(:,:,:,:,kc,jc) 
     &                    + blli(:,:,:,:,ic,kc) * tmp2r(:,:,:,:,kc,jc)
            
                  end do
               end do
            end do
            bllr = tmp1r(:,:,:,:,:,:)
            blli = tmp1i(:,:,:,:,:,:)
            tmp1r = 0.0d0
            tmp1i = 0.0d0
            do ic = 1,nc
               do jc = 1,nc
                  do kc = 1,nc
               
                     tmp1r(:,:,:,:,ic,jc) = tmp1r(:,:,:,:,ic,jc) 
     &                    + brlr(:,:,:,:,ic,kc) * tmp3r(:,:,:,:,kc,jc) 
     &                    - brli(:,:,:,:,ic,kc) * tmp3i(:,:,:,:,kc,jc)
                      
                     tmp1i(:,:,:,:,ic,jc) = tmp1i(:,:,:,:,ic,jc)  
     &                    + brlr(:,:,:,:,ic,kc) * tmp3i(:,:,:,:,kc,jc) 
     &                    + brli(:,:,:,:,ic,kc) * tmp3r(:,:,:,:,kc,jc)
                        
                  end do
               end do
            end do
            brlr = tmp1r(:,:,:,:,:,:)
            brli = tmp1i(:,:,:,:,:,:)

            if(ib < 16) then
c
c     --------------------------------------------
c     CASE 2 - quarks at (0,3), (-2,-2) and (2,-2)
c     --------------------------------------------
c     
c     top leg
c     
               tmp1r = 0.0d0
               tmp1i = 0.0d0
               sxlr = cshift(upyr(:,:,:,:,:,:),dim=yhat,shift=1)
               sxli = cshift(upyi(:,:,:,:,:,:),dim=yhat,shift=1)
         
               do ic = 1,nc
                  do jc = 1,nc
                     do kc = 1,nc
               
                        tmp1r(:,:,:,:,ic,jc) = tmp1r(:,:,:,:,ic,jc) 
     &                       + bxlr(:,:,:,:,ic,kc) * sxlr(:,:,:,:,kc,jc) - bxli(:,:,:,:,ic,kc) * sxli(:,:,:,:,kc,jc)
                
                        tmp1i(:,:,:,:,ic,jc) = tmp1i(:,:,:,:,ic,jc)  
     &                       + bxlr(:,:,:,:,ic,kc) * sxli(:,:,:,:,kc,jc) + bxli(:,:,:,:,ic,kc) * sxlr(:,:,:,:,kc,jc)
             
                     end do
                  end do
               end do
         
               bxlr = tmp1r(:,:,:,:,:,:)
               bxli = tmp1i(:,:,:,:,:,:)
               tmp1r = 0.0d0
               tmp1i = 0.0d0
               ustr = cshift(uptr(:,:,:,:,:,:),dim=yhat,shift=3) 
               usti = cshift(upti(:,:,:,:,:,:),dim=yhat,shift=3)
         
               include 'Yfiles/topleg2.f'
c
c     left leg, i.e., quark at (-2,-2)
c
               tmp1r = 0.0d0
               tmp1i = 0.0d0
               ustr = cshift(cshift(uptr(:,:,:,:,:,:),dim=xhat,shift=-2),dim=yhat,shift=-2) 
               usti = cshift(cshift(upti(:,:,:,:,:,:),dim=xhat,shift=-2),dim=yhat,shift=-2)
         
               include 'Yfiles/leftleg2.f'
c
c     right leg, i.e., quark at (2,-2)
c
               tmp1r = 0.0d0
               tmp1i = 0.0d0
               ustr = cshift(cshift(uptr(:,:,:,:,:,:),dim=xhat,shift=2),dim=yhat,shift=-2) 
               usti = cshift(cshift(upti(:,:,:,:,:,:),dim=xhat,shift=2),dim=yhat,shift=-2)
         
               include 'Yfiles/rightleg2.f'
            
               include 'Yfiles/res.f'

               W(:,:,:,:,2,it) = W(:,:,:,:,2,it) + Res(:,:,:,:) / 6.0d0
            
            endif

            tmp2r = cshift(cshift(ld2r(:,:,:,:,BD(ib,3),:,:),dim=xhat,shift=-2),dim=yhat,shift=-2)
            tmp2i = cshift(cshift(ld2i(:,:,:,:,BD(ib,3),:,:),dim=xhat,shift=-2),dim=yhat,shift=-2)
            tmp3r = cshift(cshift(rd2r(:,:,:,:,BD(ib,7),:,:),dim=xhat,shift= 2),dim=yhat,shift=-2)
            tmp3i = cshift(cshift(rd2i(:,:,:,:,BD(ib,7),:,:),dim=xhat,shift= 2),dim=yhat,shift=-2)
            tmp1r = 0.0d0
            tmp1i = 0.0d0
            do ic = 1,nc
               do jc = 1,nc
                  do kc = 1,nc
                     
                     tmp1r(:,:,:,:,ic,jc) = tmp1r(:,:,:,:,ic,jc) 
     &                    + bllr(:,:,:,:,ic,kc) * tmp2r(:,:,:,:,kc,jc) 
     &                    - blli(:,:,:,:,ic,kc) * tmp2i(:,:,:,:,kc,jc)
                
                     tmp1i(:,:,:,:,ic,jc) = tmp1i(:,:,:,:,ic,jc)
     &                    + bllr(:,:,:,:,ic,kc) * tmp2i(:,:,:,:,kc,jc) 
     &                    + blli(:,:,:,:,ic,kc) * tmp2r(:,:,:,:,kc,jc)
            
                  end do
               end do
            end do
            bllr = tmp1r(:,:,:,:,:,:)
            blli = tmp1i(:,:,:,:,:,:)
            tmp1r = 0.0d0
            tmp1i = 0.0d0
            do ic = 1,nc
               do jc = 1,nc
                  do kc = 1,nc
               
                     tmp1r(:,:,:,:,ic,jc) = tmp1r(:,:,:,:,ic,jc) 
     &                    + brlr(:,:,:,:,ic,kc) * tmp3r(:,:,:,:,kc,jc) 
     &                    - brli(:,:,:,:,ic,kc) * tmp3i(:,:,:,:,kc,jc)
                      
                     tmp1i(:,:,:,:,ic,jc) = tmp1i(:,:,:,:,ic,jc)
     &                    + brlr(:,:,:,:,ic,kc) * tmp3i(:,:,:,:,kc,jc) 
     &                    + brli(:,:,:,:,ic,kc) * tmp3r(:,:,:,:,kc,jc)
                        
                  end do
               end do
            end do
            brlr = tmp1r(:,:,:,:,:,:)
            brli = tmp1i(:,:,:,:,:,:)

            if(ib < 64) then
c
c     --------------------------------------------
c     CASE 3 - quarks at (0,5), (-4,-3) and (4,-3)
c     --------------------------------------------
c     
c     top leg
c     
               tmp1r = 0.0d0
               tmp1i = 0.0d0
               sxlr = cshift(upyr(:,:,:,:,:,:),dim=yhat,shift=3)
               sxli = cshift(upyi(:,:,:,:,:,:),dim=yhat,shift=3)
         
               do ic = 1,nc
                  do jc = 1,nc
                     do kc = 1,nc
               
                        tmp1r(:,:,:,:,ic,jc) = tmp1r(:,:,:,:,ic,jc) 
     &                       + bxlr(:,:,:,:,ic,kc) * sxlr(:,:,:,:,kc,jc) - bxli(:,:,:,:,ic,kc) * sxli(:,:,:,:,kc,jc)
           
                        tmp1i(:,:,:,:,ic,jc) = tmp1i(:,:,:,:,ic,jc)  
     &                       + bxlr(:,:,:,:,ic,kc) * sxli(:,:,:,:,kc,jc) + bxli(:,:,:,:,ic,kc) * sxlr(:,:,:,:,kc,jc)
             
                     end do
                  end do
               end do
         
               bxlr = tmp1r(:,:,:,:,:,:)
               bxli = tmp1i(:,:,:,:,:,:)
               tmp1r = 0.0d0
               tmp1i = 0.0d0
               ustr = cshift(uptr(:,:,:,:,:,:),dim=yhat,shift=5) 
               usti = cshift(upti(:,:,:,:,:,:),dim=yhat,shift=5)
         
               include 'Yfiles/topleg2.f'
c
c     left leg, i.e., quark at (-4,-3)
c
               tmp1r = 0.0d0
               tmp1i = 0.0d0

               ustr = cshift(cshift(uptr(:,:,:,:,:,:),dim=xhat,shift=-4),dim=yhat,shift=-3) 
               usti = cshift(cshift(upti(:,:,:,:,:,:),dim=xhat,shift=-4),dim=yhat,shift=-3)
         
               include 'Yfiles/leftleg2.f'
c     
c     right leg, i.e., quark at (4,-3)
c
               tmp1r = 0.0d0
               tmp1i = 0.0d0
               ustr = cshift(cshift(uptr(:,:,:,:,:,:),dim=xhat,shift=4),dim=yhat,shift=-3) 
               usti = cshift(cshift(upti(:,:,:,:,:,:),dim=xhat,shift=4),dim=yhat,shift=-3)
         
               include 'Yfiles/rightleg2.f'

               include 'Yfiles/res.f'         

               W(:,:,:,:,3,it) = W(:,:,:,:,3,it) + Res(:,:,:,:) / 6.0d0

            endif
c
c     --------------------------------------------
c     CASE 4 - quarks at (0,7), (-6,-4) and (6,-4)
c     --------------------------------------------
c     
c     top leg
c
            tmp1r = 0.0d0
            tmp1i = 0.0d0
            sxlr = cshift(upyr(:,:,:,:,:,:),dim=yhat,shift=5)
            sxli = cshift(upyi(:,:,:,:,:,:),dim=yhat,shift=5)
         
            do ic = 1,nc
               do jc = 1,nc
                  do kc = 1,nc

                     tmp1r(:,:,:,:,ic,jc) = tmp1r(:,:,:,:,ic,jc) 
     &                    + bxlr(:,:,:,:,ic,kc) * sxlr(:,:,:,:,kc,jc) - bxli(:,:,:,:,ic,kc) * sxli(:,:,:,:,kc,jc)
                  
                     tmp1i(:,:,:,:,ic,jc) = tmp1i(:,:,:,:,ic,jc)  
     &                    + bxlr(:,:,:,:,ic,kc) * sxli(:,:,:,:,kc,jc) + bxli(:,:,:,:,ic,kc) * sxlr(:,:,:,:,kc,jc)
          
                  end do
               end do
            end do
         
            bxlr = tmp1r(:,:,:,:,:,:)
            bxli = tmp1i(:,:,:,:,:,:)
            tmp1r = 0.0d0
            tmp1i = 0.0d0
            ustr = cshift(uptr(:,:,:,:,:,:),dim=yhat,shift=7) 
            usti = cshift(upti(:,:,:,:,:,:),dim=yhat,shift=7)
         
            include 'Yfiles/topleg2.f'

c     
c     left leg, i.e., quark at (-6,-4)
c
            tmp1r = 0.0d0
            tmp1i = 0.0d0
            tmp2r = cshift(cshift(ld2r(:,:,:,:,BD(ib,4),:,:),dim=xhat,shift=-4),dim=yhat,shift=-3)
            tmp2i = cshift(cshift(ld2i(:,:,:,:,BD(ib,4),:,:),dim=xhat,shift=-4),dim=yhat,shift=-3)
         
            do ic = 1,nc
               do jc = 1,nc
                  do kc = 1,nc
                  
                     tmp1r(:,:,:,:,ic,jc) = tmp1r(:,:,:,:,ic,jc) 
     &                    + bllr(:,:,:,:,ic,kc) * tmp2r(:,:,:,:,kc,jc) 
     &                    - blli(:,:,:,:,ic,kc) * tmp2i(:,:,:,:,kc,jc)

                     tmp1i(:,:,:,:,ic,jc) = tmp1i(:,:,:,:,ic,jc)
     &                    + bllr(:,:,:,:,ic,kc) * tmp2i(:,:,:,:,kc,jc) 
     &                    + blli(:,:,:,:,ic,kc) * tmp2r(:,:,:,:,kc,jc)

                  end do
               end do
            end do
         
            bllr = tmp1r(:,:,:,:,:,:)
            blli = tmp1i(:,:,:,:,:,:)
            tmp1r = 0.0d0
            tmp1i = 0.0d0
            ustr = cshift(cshift(uptr(:,:,:,:,:,:),dim=xhat,shift=-6),dim=yhat,shift=-4) 
            usti = cshift(cshift(upti(:,:,:,:,:,:),dim=xhat,shift=-6),dim=yhat,shift=-4)
         
            include 'Yfiles/leftleg2.f'
c
c     right leg, i.e., quark at (6,-4)
c
            tmp1r = 0.0d0
            tmp1i = 0.0d0
            tmp3r = cshift(cshift(rd2r(:,:,:,:,BD(ib,8),:,:),dim=xhat,shift=4),dim=yhat,shift=-3)
            tmp3i = cshift(cshift(rd2i(:,:,:,:,BD(ib,8),:,:),dim=xhat,shift=4),dim=yhat,shift=-3)
         
            do ic = 1,nc
               do jc = 1,nc
                  do kc = 1,nc
               
                     tmp1r(:,:,:,:,ic,jc) = tmp1r(:,:,:,:,ic,jc) 
     &                    + brlr(:,:,:,:,ic,kc) * tmp3r(:,:,:,:,kc,jc) 
     &                    - brli(:,:,:,:,ic,kc) * tmp3i(:,:,:,:,kc,jc)
                  
                     tmp1i(:,:,:,:,ic,jc) = tmp1i(:,:,:,:,ic,jc)
     &                    + brlr(:,:,:,:,ic,kc) * tmp3i(:,:,:,:,kc,jc) 
     &                    + brli(:,:,:,:,ic,kc) * tmp3r(:,:,:,:,kc,jc)

                  end do
               end do
            end do
         
            brlr = tmp1r(:,:,:,:,:,:)
            brli = tmp1i(:,:,:,:,:,:)
            tmp1r = 0.0d0
            tmp1i = 0.0d0
            ustr = cshift(cshift(uptr(:,:,:,:,:,:),dim=xhat,shift=6),dim=yhat,shift=-4) 
            usti = cshift(cshift(upti(:,:,:,:,:,:),dim=xhat,shift=6),dim=yhat,shift=-4)
         
            include 'Yfiles/rightleg2.f'
            
            include 'Yfiles/res.f'         

            W(:,:,:,:,4,it) = Res(:,:,:,:) / 6.0d0
      
         end do
      end do
      
      return
      
      end subroutine y_loops

      END MODULE YLOOPSSELECT
