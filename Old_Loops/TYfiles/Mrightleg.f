         do ic = 1,nc
            do jc = 1,nc
               do kc = 1,nc
               
                  tmpr(:,:,:,:,ic,jc) = tmpr(:,:,:,:,ic,jc) 
     &                 + mbrlr(:,:,:,:,kc,ic) * ustr(:,:,:,:,kc,jc) + mbrli(:,:,:,:,kc,ic) * usti(:,:,:,:,kc,jc)
              
                  tmpi(:,:,:,:,ic,jc) = tmpi(:,:,:,:,ic,jc)  
     &                 + mbrlr(:,:,:,:,kc,ic) * usti(:,:,:,:,kc,jc) - mbrli(:,:,:,:,kc,ic) * ustr(:,:,:,:,kc,jc)
        
               end do
            end do
         end do
         
         srlr = 0.0d0
         srli = 0.0d0
         tlr = cshift(mbrlr(:,:,:,:,:,:),dim=that,shift=it)
         tli = cshift(mbrli(:,:,:,:,:,:),dim=that,shift=it)
         
         do ic = 1,nc
            do jc = 1,nc
               do kc = 1,nc
               
                  srlr(:,:,:,:,ic,jc) = srlr(:,:,:,:,ic,jc) 
     &                 + tmpr(:,:,:,:,ic,kc) * tlr(:,:,:,:,kc,jc) - tmpi(:,:,:,:,ic,kc) * tli(:,:,:,:,kc,jc)
            
                  srli(:,:,:,:,ic,jc) = srli(:,:,:,:,ic,jc) 
     &                 + tmpr(:,:,:,:,ic,kc) * tli(:,:,:,:,kc,jc) + tmpi(:,:,:,:,ic,kc) * tlr(:,:,:,:,kc,jc)
           
               end do
            end do
         end do
