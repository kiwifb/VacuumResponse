         tmp1r = 0.0d0
         tmp1i = 0.0d0

         do ic = 1,nc
            do jc = 1,nc
               do kc = 1,nc

                  tmp1r(:,:,:,:,ic,jc) = tmp1r(:,:,:,:,ic,jc)
     &                 + bllr(:,:,:,:,kc,ic) * t_tr(:,:,:,:,kc,jc) + blli(:,:,:,:,kc,ic) * t_ti(:,:,:,:,kc,jc)

                  tmp1i(:,:,:,:,ic,jc) = tmp1i(:,:,:,:,ic,jc)
     &                 + bllr(:,:,:,:,kc,ic) * t_ti(:,:,:,:,kc,jc) - blli(:,:,:,:,kc,ic) * t_tr(:,:,:,:,kc,jc)

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
     &                 + tmp1r(:,:,:,:,ic,kc) * sbllr(:,:,:,:,kc,jc) - tmp1i(:,:,:,:,ic,kc) * sblli(:,:,:,:,kc,jc)

                  slli(:,:,:,:,ic,jc) = slli(:,:,:,:,ic,jc)
     &                 + tmp1r(:,:,:,:,ic,kc) * sblli(:,:,:,:,kc,jc) + tmp1i(:,:,:,:,ic,kc) * sbllr(:,:,:,:,kc,jc)

               end do
            end do
         end do
