         tmpr = 0.0d0
         tmpi = 0.0d0

         do ic = 1,nc
            do jc = 1,nc
               do kc = 1,nc

                  tmpr(:,:,:,:,ic,jc) = tmpr(:,:,:,:,ic,jc)
     &                 + bllr(:,:,:,:,ic,kc) * ustr(:,:,:,:,kc,jc) - blli(:,:,:,:,ic,kc) * usti(:,:,:,:,kc,jc)

                  tmpi(:,:,:,:,ic,jc) = tmpi(:,:,:,:,ic,jc)
     &                 + bllr(:,:,:,:,ic,kc) * usti(:,:,:,:,kc,jc) + blli(:,:,:,:,ic,kc) * ustr(:,:,:,:,kc,jc)

               end do
            end do
         end do

         sllr = 0.0d0
         slli = 0.0d0
         tlr = cshift(bllr,dim=that,shift=it)
         tli = cshift(blli,dim=that,shift=it)

         do ic = 1,nc
            do jc = 1,nc
               do kc = 1,nc

                  sllr(:,:,:,:,ic,jc) = sllr(:,:,:,:,ic,jc)
     &                 + tmpr(:,:,:,:,ic,kc) * tlr(:,:,:,:,jc,kc) + tmpi(:,:,:,:,ic,kc) * tli(:,:,:,:,jc,kc)

                  slli(:,:,:,:,ic,jc) = slli(:,:,:,:,ic,jc)
     &                 - tmpr(:,:,:,:,ic,kc) * tli(:,:,:,:,jc,kc) + tmpi(:,:,:,:,ic,kc) * tlr(:,:,:,:,jc,kc)

               end do
            end do
         end do
