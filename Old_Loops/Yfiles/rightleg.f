         tmpr = 0.0d0
         tmpi = 0.0d0

         do ic = 1,nc
            do jc = 1,nc
               do kc = 1,nc
               
                  tmpr(:,:,:,:,ic,jc) = tmpr(:,:,:,:,ic,jc) 
     &                 + brlr(:,:,:,:,ic,kc) * ustr(:,:,:,:,kc,jc) - brli(:,:,:,:,ic,kc) * usti(:,:,:,:,kc,jc)
              
                  tmpi(:,:,:,:,ic,jc) = tmpi(:,:,:,:,ic,jc)  
     &                 + brlr(:,:,:,:,ic,kc) * usti(:,:,:,:,kc,jc) + brli(:,:,:,:,ic,kc) * ustr(:,:,:,:,kc,jc)
        
               end do
            end do
         end do
         
         srlr = 0.0d0
         srli = 0.0d0
         tlr = cshift(brlr,dim=that,shift=it)
         tli = cshift(brli,dim=that,shift=it)
         
         do ic = 1,nc
            do jc = 1,nc
               do kc = 1,nc
               
                  srlr(:,:,:,:,ic,jc) = srlr(:,:,:,:,ic,jc) 
     &                 + tmpr(:,:,:,:,ic,kc) * tlr(:,:,:,:,jc,kc) + tmpi(:,:,:,:,ic,kc) * tli(:,:,:,:,jc,kc)
            
                  srli(:,:,:,:,ic,jc) = srli(:,:,:,:,ic,jc) 
     &                 - tmpr(:,:,:,:,ic,kc) * tli(:,:,:,:,jc,kc) + tmpi(:,:,:,:,ic,kc) * tlr(:,:,:,:,jc,kc)
           
               end do
            end do
         end do
