         tmpr = 0.0d0
         tmpi = 0.0d0

         do ic = 1,nc
            do jc = 1,nc
               do kc = 1,nc

                  tmpr(:,:,:,:,ic,jc) = tmpr(:,:,:,:,ic,jc)
     &                 + brlr(:,:,:,:,ic,kc) * tlr(:,:,:,:,kc,jc) - brli(:,:,:,:,ic,kc) * tli(:,:,:,:,kc,jc)

                  tmpi(:,:,:,:,ic,jc) = tmpi(:,:,:,:,ic,jc)
     &                 + brlr(:,:,:,:,ic,kc) * tli(:,:,:,:,kc,jc) + brli(:,:,:,:,ic,kc) * tlr(:,:,:,:,kc,jc)

               end do
            end do
         end do

         srlr = 0.0d0
         srli = 0.0d0
         sbrlr = cshift(brlr,dim=that,shift=it)
         sbrli = cshift(brli,dim=that,shift=it)

         do ic = 1,nc
            do jc = 1,nc
               do kc = 1,nc

                  srlr(:,:,:,:,ic,jc) = srlr(:,:,:,:,ic,jc)
     &                 + tmpr(:,:,:,:,ic,kc) * sbrlr(:,:,:,:,jc,kc) + tmpi(:,:,:,:,ic,kc) * sbrli(:,:,:,:,jc,kc)

                  srli(:,:,:,:,ic,jc) = srli(:,:,:,:,ic,jc)
     &                 - tmpr(:,:,:,:,ic,kc) * sbrli(:,:,:,:,jc,kc) + tmpi(:,:,:,:,ic,kc) * sbrlr(:,:,:,:,jc,kc)

               end do
            end do
         end do
