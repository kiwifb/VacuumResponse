         do ic = 1,nc
            do jc = 1,nc
               do kc = 1,nc
               
                  tmpr(:,:,:,:,ic,jc) = tmpr(:,:,:,:,ic,jc) 
     &                 + mbllr(:,:,:,:,kc,ic) * ustr(:,:,:,:,kc,jc) + mblli(:,:,:,:,kc,ic) * usti(:,:,:,:,kc,jc)
          
                  tmpi(:,:,:,:,ic,jc) = tmpi(:,:,:,:,ic,jc)  
     &                 + mbllr(:,:,:,:,kc,ic) * usti(:,:,:,:,kc,jc) - mblli(:,:,:,:,kc,ic) * ustr(:,:,:,:,kc,jc)
             
               end do
            end do
         end do
         
         sllr = 0.0d0
         slli = 0.0d0
         tlr = cshift(mbllr(:,:,:,:,:,:),dim=that,shift=it)
         tli = cshift(mblli(:,:,:,:,:,:),dim=that,shift=it)
         
         do ic = 1,nc
            do jc = 1,nc
               do kc = 1,nc
               
                  sllr(:,:,:,:,ic,jc) = sllr(:,:,:,:,ic,jc) 
     &                 + tmpr(:,:,:,:,ic,kc) * tlr(:,:,:,:,kc,jc) - tmpi(:,:,:,:,ic,kc) * tli(:,:,:,:,kc,jc)
      
                  slli(:,:,:,:,ic,jc) = slli(:,:,:,:,ic,jc) 
     &                 + tmpr(:,:,:,:,ic,kc) * tli(:,:,:,:,kc,jc) + tmpi(:,:,:,:,ic,kc) * tlr(:,:,:,:,kc,jc)
       
               end do
            end do
         end do
