         do ic = 1,nc
            do jc = 1,nc
               do kc = 1,nc
               
                  tmp1r(:,:,:,:,ic,jc) = tmp1r(:,:,:,:,ic,jc) 
     &                 + brlr(:,:,:,:,ic,kc) * ustr(:,:,:,:,kc,jc) - brli(:,:,:,:,ic,kc) * usti(:,:,:,:,kc,jc)
              
                  tmp1i(:,:,:,:,ic,jc) = tmp1i(:,:,:,:,ic,jc)  
     &                 + brlr(:,:,:,:,ic,kc) * usti(:,:,:,:,kc,jc) + brli(:,:,:,:,ic,kc) * ustr(:,:,:,:,kc,jc)
        
               end do
            end do
         end do
         
         srlr = 0.0d0
         srli = 0.0d0
         tlr = cshift(brlr(:,:,:,:,:,:),dim=that,shift=it)
         tli = cshift(brli(:,:,:,:,:,:),dim=that,shift=it)
         
         do ic = 1,nc
            do jc = 1,nc
               do kc = 1,nc
               
                  srlr(:,:,:,:,ic,jc) = srlr(:,:,:,:,ic,jc) 
     &                 + tmp1r(:,:,:,:,ic,kc) * tlr(:,:,:,:,jc,kc) + tmp1i(:,:,:,:,ic,kc) * tli(:,:,:,:,jc,kc)
            
                  srli(:,:,:,:,ic,jc) = srli(:,:,:,:,ic,jc) 
     &                 - tmp1r(:,:,:,:,ic,kc) * tli(:,:,:,:,jc,kc) + tmp1i(:,:,:,:,ic,kc) * tlr(:,:,:,:,jc,kc)
           
               end do
            end do
         end do
