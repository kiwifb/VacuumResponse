!
!  box-mm third quadrant
!
      do ic = 1, nc
         do jc = 1, nc
            do kc = 1, nc

               boxr(:,:,:,:,ic,jc) = boxr(:,:,:,:,ic,jc) +
     &              ux1r(:,:,:,:,kc,ic)*uy2r(:,:,:,:,jc,kc) - ux1i(:,:,:,:,kc,ic)*uy2i(:,:,:,:,jc,kc) +
     &              uy1r(:,:,:,:,kc,ic)*ux2r(:,:,:,:,jc,kc) - uy1i(:,:,:,:,kc,ic)*ux2i(:,:,:,:,jc,kc)

               boxi(:,:,:,:,ic,jc) = boxi(:,:,:,:,ic,jc) -
     &              ux1r(:,:,:,:,kc,ic)*uy2i(:,:,:,:,jc,kc) - ux1i(:,:,:,:,kc,ic)*uy2r(:,:,:,:,jc,kc) -
     &              uy1r(:,:,:,:,kc,ic)*ux2i(:,:,:,:,jc,kc) - uy1i(:,:,:,:,kc,ic)*ux2r(:,:,:,:,jc,kc)
            end do
         end do
      end do
