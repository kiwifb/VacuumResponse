c
c     ----------------------------------------------
c     CASE 5 - quarks at ( 4,0), (-3,-4) and (-3, 4)
c     ----------------------------------------------
c     
c     top leg ( 4,0)
c     
         call product(ur,ui,xdir,bxlr,bxli,4)

         ustr = cshift(uptr,dim=xdir,shift=4)
         usti = cshift(upti,dim=xdir,shift=4)
         
         include 'Yfiles/topleg.f'
c
c     left leg, i.e., quark at (-3, 4)
c
c     to (-1, 1)
c
         bllr = ld1r
         blli = ld1i
c
c     to (-2, 3)
c
         tmpr = 0.d0
         tmpi = 0.d0
         stmpr = cshift(cshift(ld2r,dim=xdir,shift=-1),dim=ydir,shift= 1)
         stmpi = cshift(cshift(ld2i,dim=xdir,shift=-1),dim=ydir,shift= 1)
         
         do ic = 1,nc
            do jc = 1,nc
               do kc = 1,nc
               
                  tmpr(:,:,:,:,ic,jc) = tmpr(:,:,:,:,ic,jc) 
     &                 + bllr(:,:,:,:,ic,kc) * stmpr(:,:,:,:,kc,jc) - blli(:,:,:,:,ic,kc) * stmpi(:,:,:,:,kc,jc)
                
                  tmpi(:,:,:,:,ic,jc) = tmpi(:,:,:,:,ic,jc)  
     &                 + bllr(:,:,:,:,ic,kc) * stmpi(:,:,:,:,kc,jc) + blli(:,:,:,:,ic,kc) * stmpr(:,:,:,:,kc,jc)
             
               end do
            end do
         end do
c
c     to (-3, 4)
c
         bllr = 0.d0
         blli = 0.d0
         stmpr = cshift(cshift(ld1r,dim=xdir,shift=-2),dim=ydir,shift= 3)
         stmpi = cshift(cshift(ld1i,dim=xdir,shift=-2),dim=ydir,shift= 3)
         
         do ic = 1,nc
            do jc = 1,nc
               do kc = 1,nc
               
                  bllr(:,:,:,:,ic,jc) = bllr(:,:,:,:,ic,jc) 
     &                 + tmpr(:,:,:,:,ic,kc) * stmpr(:,:,:,:,kc,jc) - tmpi(:,:,:,:,ic,kc) * stmpi(:,:,:,:,kc,jc)
                
                  blli(:,:,:,:,ic,jc) = blli(:,:,:,:,ic,jc)  
     &                 + tmpr(:,:,:,:,ic,kc) * stmpi(:,:,:,:,kc,jc) + tmpi(:,:,:,:,ic,kc) * stmpr(:,:,:,:,kc,jc)
             
               end do
            end do
         end do
         
         ustr = cshift(cshift(uptr,dim=xdir,shift=-3),dim=ydir,shift= 4) 
         usti = cshift(cshift(upti,dim=xdir,shift=-3),dim=ydir,shift= 4)
         
         include 'Yfiles/leftleg.f'
c     
c     right leg, i.e., quark at (-3,-4)
c
c     to (-1,-1)
c
         brlr = rd1r
         brli = rd1i
c
c     to (-2,-3)
c
         tmpr = 0.d0
         tmpi = 0.d0
         stmpr = cshift(cshift(rd2r,dim=xdir,shift=-1),dim=ydir,shift=-1)
         stmpi = cshift(cshift(rd2i,dim=xdir,shift=-1),dim=ydir,shift=-1)
         
         do ic = 1,nc
            do jc = 1,nc
               do kc = 1,nc
               
                  tmpr(:,:,:,:,ic,jc) = tmpr(:,:,:,:,ic,jc) 
     &                 + brlr(:,:,:,:,ic,kc) * stmpr(:,:,:,:,kc,jc) - brli(:,:,:,:,ic,kc) * stmpi(:,:,:,:,kc,jc)
                
                  tmpi(:,:,:,:,ic,jc) = tmpi(:,:,:,:,ic,jc)  
     &                 + brlr(:,:,:,:,ic,kc) * stmpi(:,:,:,:,kc,jc) + brli(:,:,:,:,ic,kc) * stmpr(:,:,:,:,kc,jc)
             
               end do
            end do
         end do
c
c     to (-3,-4)
c
         brlr = 0.d0
         brli = 0.d0
         stmpr = cshift(cshift(rd1r,dim=xdir,shift=-2),dim=ydir,shift=-3)
         stmpi = cshift(cshift(rd1i,dim=xdir,shift=-2),dim=ydir,shift=-3)
         
         do ic = 1,nc
            do jc = 1,nc
               do kc = 1,nc
               
                  brlr(:,:,:,:,ic,jc) = brlr(:,:,:,:,ic,jc) 
     &                 + tmpr(:,:,:,:,ic,kc) * stmpr(:,:,:,:,kc,jc) - tmpi(:,:,:,:,ic,kc) * stmpi(:,:,:,:,kc,jc)
                
                  brli(:,:,:,:,ic,jc) = brli(:,:,:,:,ic,jc)  
     &                 + tmpr(:,:,:,:,ic,kc) * stmpi(:,:,:,:,kc,jc) + tmpi(:,:,:,:,ic,kc) * stmpr(:,:,:,:,kc,jc)
             
               end do
            end do
         end do

         ustr = cshift(cshift(uptr,dim=xdir,shift=-3),dim=ydir,shift=-4) 
         usti = cshift(cshift(upti,dim=xdir,shift=-3),dim=ydir,shift=-4)
         
         include 'Yfiles/rightleg.f'

         include 'Yfiles/res.f'
