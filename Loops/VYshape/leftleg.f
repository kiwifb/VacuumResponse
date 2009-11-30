         tmp1r = 0.0d0
         tmp1i = 0.0d0

         do ic = 1,nc
            do jc = 1,nc
               do kc = 1,nc

                  tmp1r(:,:,:,:,ic,jc) = tmp1r(:,:,:,:,ic,jc)
     &                 + bllr(:,:,:,:,ic,kc) * t_ur(:,:,:,:,kc,jc) - blli(:,:,:,:,ic,kc) * t_ui(:,:,:,:,kc,jc)

                  tmp1i(:,:,:,:,ic,jc) = tmp1i(:,:,:,:,ic,jc)
     &                 + bllr(:,:,:,:,ic,kc) * t_ui(:,:,:,:,kc,jc) + blli(:,:,:,:,ic,kc) * t_ur(:,:,:,:,kc,jc)

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
     &                 + tmp1r(:,:,:,:,ic,kc) * sbllr(:,:,:,:,jc,kc) + tmp1i(:,:,:,:,ic,kc) * sblli(:,:,:,:,jc,kc)

                  slli(:,:,:,:,ic,jc) = slli(:,:,:,:,ic,jc)
     &                 - tmp1r(:,:,:,:,ic,kc) * sblli(:,:,:,:,jc,kc) + tmp1i(:,:,:,:,ic,kc) * sbllr(:,:,:,:,jc,kc)

               end do
            end do
         end do
