         do ic = 1,nc
            do jc = 1,nc
               do kc = 1,nc
               
                  tmpr(:,:,:,:,ic,jc) = tmpr(:,:,:,:,ic,jc) 
     &                 + mbxlr(:,:,:,:,kc,ic) * ustr(:,:,:,:,kc,jc) + mbxli(:,:,:,:,kc,ic) * usti(:,:,:,:,kc,jc)
         
                  tmpi(:,:,:,:,ic,jc) = tmpi(:,:,:,:,ic,jc)  
     &                 + mbxlr(:,:,:,:,kc,ic) * usti(:,:,:,:,kc,jc) - mbxli(:,:,:,:,kc,ic) * ustr(:,:,:,:,kc,jc)
            
               end do
            end do
         end do

         sxlr = 0.0d0
         sxli = 0.0d0
         tlr = cshift(mbxlr(:,:,:,:,:,:),dim=that,shift=it) 
         tli = cshift(mbxli(:,:,:,:,:,:),dim=that,shift=it) 
         
         do ic = 1,nc
            do jc = 1,nc
               do kc = 1,nc

                  sxlr(:,:,:,:,ic,jc) = sxlr(:,:,:,:,ic,jc) 
     &                 + tmpr(:,:,:,:,ic,kc) * tlr(:,:,:,:,kc,jc) - tmpi(:,:,:,:,ic,kc) * tli(:,:,:,:,kc,jc)
            
                  sxli(:,:,:,:,ic,jc) = sxli(:,:,:,:,ic,jc) 
     &                 + tmpr(:,:,:,:,ic,kc) * tli(:,:,:,:,kc,jc) + tmpi(:,:,:,:,ic,kc) * tlr(:,:,:,:,kc,jc)
          
               end do
            end do
         end do
