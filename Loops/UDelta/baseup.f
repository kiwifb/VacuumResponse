      bllr = 0.0d0
      blli = 0.0d0

      do ic = 1,nc
         do jc = 1,nc
            do kc = 1,nc

               bllr(:,:,:,:,ic,jc) = bllr(:,:,:,:,ic,jc)
     &              + tmp1r(:,:,:,:,ic,kc) * dylr(:,:,:,:,jc,kc) + tmp1i(:,:,:,:,ic,kc) * dyli(:,:,:,:,jc,kc)

               blli(:,:,:,:,ic,jc) = blli(:,:,:,:,ic,jc)
     &              - tmp1r(:,:,:,:,ic,kc) * dyli(:,:,:,:,jc,kc) + tmp1i(:,:,:,:,ic,kc) * dylr(:,:,:,:,jc,kc)

            end do
         end do
      end do
