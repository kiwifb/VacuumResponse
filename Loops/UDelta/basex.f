      bxlr = 0.0d0
      bxli = 0.0d0

      do ic = 1,nc
         do jc = 1,nc
            do kc = 1,nc

               bxlr(:,:,:,:,ic,jc) = bxlr(:,:,:,:,ic,jc)
     &              + dxlr(:,:,:,:,ic,kc) * cxlr(:,:,:,:,kc,jc) - dxli(:,:,:,:,ic,kc) * cxli(:,:,:,:,kc,jc)

               bxli(:,:,:,:,ic,jc) = bxli(:,:,:,:,ic,jc)
     &              + dxlr(:,:,:,:,ic,kc) * cxli(:,:,:,:,kc,jc) + dxli(:,:,:,:,ic,kc) * cxlr(:,:,:,:,kc,jc)

            end do
         end do
      end do
