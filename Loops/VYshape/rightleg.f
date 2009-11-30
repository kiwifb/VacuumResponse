         tmp1r = 0.0d0
         tmp1i = 0.0d0

         do ic = 1,nc
            do jc = 1,nc
               do kc = 1,nc

                  tmp1r(:,:,:,:,ic,jc) = tmp1r(:,:,:,:,ic,jc)
     &                 + brlr(:,:,:,:,ic,kc) * t_dr(:,:,:,:,kc,jc) - brli(:,:,:,:,ic,kc) * t_di(:,:,:,:,kc,jc)

                  tmp1i(:,:,:,:,ic,jc) = tmp1i(:,:,:,:,ic,jc)
     &                 + brlr(:,:,:,:,ic,kc) * t_di(:,:,:,:,kc,jc) + brli(:,:,:,:,ic,kc) * t_dr(:,:,:,:,kc,jc)

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
     &                 + tmp1r(:,:,:,:,ic,kc) * sbrlr(:,:,:,:,jc,kc) + tmp1i(:,:,:,:,ic,kc) * sbrli(:,:,:,:,jc,kc)

                  srli(:,:,:,:,ic,jc) = srli(:,:,:,:,ic,jc)
     &                 - tmp1r(:,:,:,:,ic,kc) * sbrli(:,:,:,:,jc,kc) + tmp1i(:,:,:,:,ic,kc) * sbrlr(:,:,:,:,jc,kc)

               end do
            end do
         end do
