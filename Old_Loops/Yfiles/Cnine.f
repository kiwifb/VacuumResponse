c
c     ----------------------------------------------
c     CASE 9 - quarks at ( 9,0), (-5,-8) and (-5, 8)
c     ----------------------------------------------
c     
c     top leg ( 9,0)
c     
         call product(ur,ui,xdir,bxlr,bxli,9)

         ustr = cshift(uptr,dim=xdir,shift=9)
         usti = cshift(upti,dim=xdir,shift=9)
         
         include 'Yfiles/topleg.f'
c
c     left leg, i.e., quark at (-5, 8)
c
c     to (-2, 3)
c
         bllr = ld3r
         blli = ld3i
c
c     to (-3, 5)
c
         tmpr = 0.d0
         tmpi = 0.d0
         stmpr = cshift(cshift(ld2r,dim=xdir,shift=-2),dim=ydir,shift= 3)
         stmpi = cshift(cshift(ld2i,dim=xdir,shift=-2),dim=ydir,shift= 3)
         
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
c     to (-5, 8)
c
         bllr = 0.d0
         blli = 0.d0
         stmpr = cshift(cshift(ld3r,dim=xdir,shift=-3),dim=ydir,shift= 5)
         stmpi = cshift(cshift(ld3i,dim=xdir,shift=-3),dim=ydir,shift= 5)
         
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

         ustr = cshift(cshift(uptr,dim=xdir,shift=-5),dim=ydir,shift= 8) 
         usti = cshift(cshift(upti,dim=xdir,shift=-5),dim=ydir,shift= 8)
         
         include 'Yfiles/leftleg.f'
c     
c     right leg, i.e., quark at (-5,-8)
c
c     to (-2,-3)
c
         brlr = rd3r
         brli = rd3i
c
c     to (-3,-5)
c
         tmpr = 0.d0
         tmpi = 0.d0
         stmpr = cshift(cshift(rd2r,dim=xdir,shift=-2),dim=ydir,shift=-3)
         stmpi = cshift(cshift(rd2i,dim=xdir,shift=-2),dim=ydir,shift=-3)
         
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
c     to (-5,-8)
c
         brlr = 0.d0
         brli = 0.d0
         stmpr = cshift(cshift(rd3r,dim=xdir,shift=-3),dim=ydir,shift=-5)
         stmpi = cshift(cshift(rd3i,dim=xdir,shift=-3),dim=ydir,shift=-5)
         
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
         
         ustr = cshift(cshift(uptr,dim=xdir,shift=-5),dim=ydir,shift=-8) 
         usti = cshift(cshift(upti,dim=xdir,shift=-5),dim=ydir,shift=-8)
         
         include 'Yfiles/rightleg.f'

         include 'Yfiles/res.f'
