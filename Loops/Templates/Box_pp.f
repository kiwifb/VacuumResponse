!
!  box-pp first quadrant
!
      do ic = 1, nc
         do jc = 1, nc
            do kc = 1, nc

               boxr(:,:,:,:,ic,jc) = boxr(:,:,:,:,ic,jc) +
     &              ux1r(:,:,:,:,ic,kc)*uy2r(:,:,:,:,kc,jc) - ux1i(:,:,:,:,ic,kc)*uy2i(:,:,:,:,kc,jc) +
     &              uy1r(:,:,:,:,ic,kc)*ux2r(:,:,:,:,kc,jc) - uy1i(:,:,:,:,ic,kc)*ux2i(:,:,:,:,kc,jc)

               boxi(:,:,:,:,ic,jc) = boxi(:,:,:,:,ic,jc) +
     &              ux1r(:,:,:,:,ic,kc)*uy2i(:,:,:,:,kc,jc) + ux1i(:,:,:,:,ic,kc)*uy2r(:,:,:,:,kc,jc) +
     &              uy1r(:,:,:,:,ic,kc)*ux2i(:,:,:,:,kc,jc) + uy1i(:,:,:,:,ic,kc)*ux2r(:,:,:,:,kc,jc)

            end do
         end do
      end do
