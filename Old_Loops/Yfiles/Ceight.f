c
c     ----------------------------------------------
c     CASE 8 - quarks at ( 8,0), (-4,-7) and (-4, 7)
c     ----------------------------------------------
c     
c     top leg ( 8,0)
c     
         call product(ur,ui,xdir,bxlr,bxli,8)

         ustr = cshift(uptr,dim=xdir,shift=8)
         usti = cshift(upti,dim=xdir,shift=8)
         
         include 'Yfiles/topleg.f'
c
c     left leg, i.e., quark at (-4, 7)
c
c     to (-1, 2)
c
         bllr = ld2r
         blli = ld2i
c
c     to (-3, 5)
c
         tmpr = 0.d0
         tmpi = 0.d0
         stmpr = cshift(cshift(ld3r,dim=xdir,shift=-1),dim=ydir,shift= 2)
         stmpi = cshift(cshift(ld3i,dim=xdir,shift=-1),dim=ydir,shift= 2)
         
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
c     to (-4, 7)
c
         bllr = 0.d0
         blli = 0.d0
         stmpr = cshift(cshift(ld2r,dim=xdir,shift=-3),dim=ydir,shift= 5)
         stmpi = cshift(cshift(ld2i,dim=xdir,shift=-3),dim=ydir,shift= 5)
         
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

         ustr = cshift(cshift(uptr,dim=xdir,shift=-4),dim=ydir,shift= 7) 
         usti = cshift(cshift(upti,dim=xdir,shift=-4),dim=ydir,shift= 7)
         
         include 'Yfiles/leftleg.f'
c     
c     right leg, i.e., quark at (-4,-7)
c
c     to (-1,-2)
c
         brlr = rd2r
         brli = rd2i
c
c     to (-3,-5)
c
         tmpr = 0.d0
         tmpi = 0.d0
         stmpr = cshift(cshift(rd3r,dim=xdir,shift=-1),dim=ydir,shift=-2)
         stmpi = cshift(cshift(rd3i,dim=xdir,shift=-1),dim=ydir,shift=-2)
         
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
c     to (-4,-7)
c
         brlr = 0.d0
         brli = 0.d0
         stmpr = cshift(cshift(rd2r,dim=xdir,shift=-3),dim=ydir,shift=-5)
         stmpi = cshift(cshift(rd2i,dim=xdir,shift=-3),dim=ydir,shift=-5)
         
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
         
         ustr = cshift(cshift(uptr,dim=xdir,shift=-4),dim=ydir,shift=-7) 
         usti = cshift(cshift(upti,dim=xdir,shift=-4),dim=ydir,shift=-7)
         
         include 'Yfiles/rightleg.f'

         include 'Yfiles/res.f'
