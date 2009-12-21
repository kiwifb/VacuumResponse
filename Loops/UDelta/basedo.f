      brlr = 0.0d0
      brli = 0.0d0

      do ic = 1,nc
         do jc = 1,nc
            do kc = 1,nc

               brlr(:,:,:,:,ic,jc) = brlr(:,:,:,:,ic,jc)
     &              + tmp1r(:,:,:,:,kc,ic) * dylr(:,:,:,:,jc,kc) - tmp1i(:,:,:,:,kc,ic) * dyli(:,:,:,:,jc,kc)

               brli(:,:,:,:,ic,jc) = brli(:,:,:,:,ic,jc)
     &              - tmp1r(:,:,:,:,kc,ic) * dyli(:,:,:,:,jc,kc) - tmp1i(:,:,:,:,kc,ic) * dylr(:,:,:,:,jc,kc)

            end do
         end do
      end do
