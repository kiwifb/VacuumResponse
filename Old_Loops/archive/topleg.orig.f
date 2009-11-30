         do ic = 1,nc
            do jc = 1,nc
               do kc = 1,nc
               
                  tmpr(:,:,:,:,ic,jc) = tmpr(:,:,:,:,ic,jc) 
     &                 + bylr(:,:,:,:,ic,kc) * ustr(:,:,:,:,kc,jc) - byli(:,:,:,:,ic,kc) * usti(:,:,:,:,kc,jc)
         
                  tmpi(:,:,:,:,ic,jc) = tmpi(:,:,:,:,ic,jc)  
     &                 + bylr(:,:,:,:,ic,kc) * usti(:,:,:,:,kc,jc) + byli(:,:,:,:,ic,kc) * ustr(:,:,:,:,kc,jc)
            
               end do
            end do
         end do

         sylr = 0.0d0
         syli = 0.0d0
         tlr = cshift(bylr(:,:,:,:,:,:),dim=that,shift=it) 
         tli = cshift(byli(:,:,:,:,:,:),dim=that,shift=it) 
         
         do ic = 1,nc
            do jc = 1,nc
               do kc = 1,nc

                  sylr(:,:,:,:,ic,jc) = sylr(:,:,:,:,ic,jc) 
     &                 + tmpr(:,:,:,:,ic,kc) * tlr(:,:,:,:,jc,kc) + tmpi(:,:,:,:,ic,kc) * tli(:,:,:,:,jc,kc)
            
                  syli(:,:,:,:,ic,jc) = syli(:,:,:,:,ic,jc) 
     &                 - tmpr(:,:,:,:,ic,kc) * tli(:,:,:,:,jc,kc) + tmpi(:,:,:,:,ic,kc) * tlr(:,:,:,:,jc,kc)
          
               end do
            end do
         end do
