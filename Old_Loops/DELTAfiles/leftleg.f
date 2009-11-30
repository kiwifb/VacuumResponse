         tmpr = 0.0d0
         tmpi = 0.0d0

         do ic = 1,nc
            do jc = 1,nc
               do kc = 1,nc

                  tmpr(:,:,:,:,ic,jc) = tmpr(:,:,:,:,ic,jc)
     &                 + bllr(:,:,:,:,ic,kc) * tl1r(:,:,:,:,kc,jc) - blli(:,:,:,:,ic,kc) * tl1i(:,:,:,:,kc,jc)

                  tmpi(:,:,:,:,ic,jc) = tmpi(:,:,:,:,ic,jc)
     &                 + bllr(:,:,:,:,ic,kc) * tl1i(:,:,:,:,kc,jc) + blli(:,:,:,:,ic,kc) * tl1r(:,:,:,:,kc,jc)

               end do
            end do
         end do

         sllr = 0.0d0
         slli = 0.0d0
         sbllr = cshift(bllr,dim=that,shift=it)
         sblli = cshift(blli,dim=that,shift=it)

         do ic = 1,nc
            do jc = 1,nc
               do kc = 1,nc

                  sllr(:,:,:,:,ic,jc) = sllr(:,:,:,:,ic,jc)
     &                 + tmpr(:,:,:,:,ic,kc) * sbllr(:,:,:,:,jc,kc) + tmpi(:,:,:,:,ic,kc) * sblli(:,:,:,:,jc,kc)

                  slli(:,:,:,:,ic,jc) = slli(:,:,:,:,ic,jc)
     &                 - tmpr(:,:,:,:,ic,kc) * sblli(:,:,:,:,jc,kc) + tmpi(:,:,:,:,ic,kc) * sbllr(:,:,:,:,jc,kc)

               end do
            end do
         end do
