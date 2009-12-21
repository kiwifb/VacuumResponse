         tmpr = 0.0d0
         tmpi = 0.0d0

         do ic = 1,nc
            do jc = 1,nc
               do kc = 1,nc

                  tmpr(:,:,:,:,ic,jc) = tmpr(:,:,:,:,ic,jc)
     &                 + bxlr(:,:,:,:,ic,kc) * ustr(:,:,:,:,kc,jc) - bxli(:,:,:,:,ic,kc) * usti(:,:,:,:,kc,jc)

                  tmpi(:,:,:,:,ic,jc) = tmpi(:,:,:,:,ic,jc)
     &                 + bxlr(:,:,:,:,ic,kc) * usti(:,:,:,:,kc,jc) + bxli(:,:,:,:,ic,kc) * ustr(:,:,:,:,kc,jc)

               end do
            end do
         end do

         sxlr = 0.0d0
         sxli = 0.0d0
         tlr = cshift(bxlr,dim=that,shift=it) 
         tli = cshift(bxli,dim=that,shift=it) 

         do ic = 1,nc
            do jc = 1,nc
               do kc = 1,nc

                  sxlr(:,:,:,:,ic,jc) = sxlr(:,:,:,:,ic,jc)
     &                 + tmpr(:,:,:,:,ic,kc) * tlr(:,:,:,:,jc,kc) + tmpi(:,:,:,:,ic,kc) * tli(:,:,:,:,jc,kc)

                  sxli(:,:,:,:,ic,jc) = sxli(:,:,:,:,ic,jc)
     &                 - tmpr(:,:,:,:,ic,kc) * tli(:,:,:,:,jc,kc) + tmpi(:,:,:,:,ic,kc) * tlr(:,:,:,:,jc,kc)

               end do
            end do
         end do
