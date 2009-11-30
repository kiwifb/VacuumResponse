!
!  box-mp second quadrant
!
      do ic = 1, nc
         do jc = 1, nc
            do kc = 1, nc

               boxr(:,:,:,:,ic,jc) = boxr(:,:,:,:,ic,jc) +
     &              ux1r(:,:,:,:,kc,ic)*uy2r(:,:,:,:,kc,jc) + ux1i(:,:,:,:,kc,ic)*uy2i(:,:,:,:,kc,jc) +
     &              uy1r(:,:,:,:,ic,kc)*ux2r(:,:,:,:,jc,kc) + uy1i(:,:,:,:,ic,kc)*ux2i(:,:,:,:,jc,kc)

               boxi(:,:,:,:,ic,jc) = boxi(:,:,:,:,ic,jc) +
     &              ux1r(:,:,:,:,kc,ic)*uy2i(:,:,:,:,kc,jc) - ux1i(:,:,:,:,kc,ic)*uy2r(:,:,:,:,kc,jc) -
     &              uy1r(:,:,:,:,ic,kc)*ux2i(:,:,:,:,jc,kc) + uy1i(:,:,:,:,ic,kc)*ux2r(:,:,:,:,jc,kc)

            end do
         end do
      end do
