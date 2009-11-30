         do ic = 1,nc
            do jc = 1,nc
               do kc = 1,nc
               
                  tmp1r(:,:,:,:,ic,jc) = tmp1r(:,:,:,:,ic,jc) 
     &                 + bllr(:,:,:,:,ic,kc) * ustr(:,:,:,:,kc,jc) - blli(:,:,:,:,ic,kc) * usti(:,:,:,:,kc,jc)
          
                  tmp1i(:,:,:,:,ic,jc) = tmp1i(:,:,:,:,ic,jc)  
     &                 + bllr(:,:,:,:,ic,kc) * usti(:,:,:,:,kc,jc) + blli(:,:,:,:,ic,kc) * ustr(:,:,:,:,kc,jc)
             
               end do
            end do
         end do
         
         sllr = 0.0d0
         slli = 0.0d0
         tlr = cshift(bllr(:,:,:,:,:,:),dim=that,shift=it)
         tli = cshift(blli(:,:,:,:,:,:),dim=that,shift=it)
         
         do ic = 1,nc
            do jc = 1,nc
               do kc = 1,nc
               
                  sllr(:,:,:,:,ic,jc) = sllr(:,:,:,:,ic,jc) 
     &                 + tmp1r(:,:,:,:,ic,kc) * tlr(:,:,:,:,jc,kc) + tmp1i(:,:,:,:,ic,kc) * tli(:,:,:,:,jc,kc)
      
                  slli(:,:,:,:,ic,jc) = slli(:,:,:,:,ic,jc) 
     &                 - tmp1r(:,:,:,:,ic,kc) * tli(:,:,:,:,jc,kc) + tmp1i(:,:,:,:,ic,kc) * tlr(:,:,:,:,jc,kc)
       
               end do
            end do
         end do
