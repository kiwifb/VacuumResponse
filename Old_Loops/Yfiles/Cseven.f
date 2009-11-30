c
c     ----------------------------------------------
c     CASE 7 - quarks at ( 7,0), (-4,-6) and (-4, 6)
c     ----------------------------------------------
c     
c     top leg ( 7,0)
c     
         call product(ur,ui,xdir,bxlr,bxli,6)
         call product(ur,ui,xdir,bxlr,bxli,7)

         ustr = cshift(uptr,dim=xdir,shift=7)
         usti = cshift(upti,dim=xdir,shift=7)
         
         include 'Yfiles/topleg.f'
c
c     left leg, i.e., quark at (-4, 6)
c
c     to (-2, 3)
c
         tmpr = ld3r
         tmpi = ld3i
c
c     to (-4, 6)
c
         bllr = 0.d0
         blli = 0.d0
         stmpr = cshift(cshift(ld3r,dim=xdir,shift=-2),dim=ydir,shift= 3)
         stmpi = cshift(cshift(ld3i,dim=xdir,shift=-2),dim=ydir,shift= 3)
         
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

         ustr = cshift(cshift(uptr,dim=xdir,shift=-4),dim=ydir,shift= 6) 
         usti = cshift(cshift(upti,dim=xdir,shift=-4),dim=ydir,shift= 6)
         
         include 'Yfiles/leftleg.f'
c     
c     right leg, i.e., quark at (-4,-6)
c
c     to (-2,-3)
c
         tmpr = rd3r
         tmpi = rd3i
c
c     to (-4,-6)
c
         brlr = 0.d0
         brli = 0.d0
         stmpr = cshift(cshift(rd3r,dim=xdir,shift=-2),dim=ydir,shift=-3)
         stmpi = cshift(cshift(rd3i,dim=xdir,shift=-2),dim=ydir,shift=-3)
         
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
         
         ustr = cshift(cshift(uptr,dim=xdir,shift=-4),dim=ydir,shift=-6) 
         usti = cshift(cshift(upti,dim=xdir,shift=-4),dim=ydir,shift=-6)
         
         include 'Yfiles/rightleg.f'

         include 'Yfiles/res.f'
