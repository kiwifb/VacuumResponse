!
!  box-pm fourth quadrant
!
      do ic = 1, nc
         do jc = 1, nc
            do kc = 1, nc

               boxr(:,:,:,:,ic,jc) = boxr(:,:,:,:,ic,jc) +
     &              ux1r(:,:,:,:,ic,kc)*uy2r(:,:,:,:,jc,kc) + ux1i(:,:,:,:,ic,kc)*uy2i(:,:,:,:,jc,kc) +
     &              uy1r(:,:,:,:,kc,ic)*ux2r(:,:,:,:,kc,jc) + uy1i(:,:,:,:,kc,ic)*ux2i(:,:,:,:,kc,jc)

               boxi(:,:,:,:,ic,jc) = boxi(:,:,:,:,ic,jc) -
     &              ux1r(:,:,:,:,ic,kc)*uy2i(:,:,:,:,jc,kc) + ux1i(:,:,:,:,ic,kc)*uy2r(:,:,:,:,jc,kc) +
     &              uy1r(:,:,:,:,kc,ic)*ux2i(:,:,:,:,kc,jc) - uy1i(:,:,:,:,kc,ic)*ux2r(:,:,:,:,kc,jc)

            end do
         end do
      end do
