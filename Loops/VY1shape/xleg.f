         tmp1r = 0.0d0
         tmp1i = 0.0d0

         do ic = 1,nc
            do jc = 1,nc
               do kc = 1,nc

                  tmp1r(:,:,:,:,ic,jc) = tmp1r(:,:,:,:,ic,jc)
     &                 + bxlr(:,:,:,:,kc,ic) * t_dr(:,:,:,:,kc,jc) + bxli(:,:,:,:,kc,ic) * t_di(:,:,:,:,kc,jc)

                  tmp1i(:,:,:,:,ic,jc) = tmp1i(:,:,:,:,ic,jc)
     &                 + bxlr(:,:,:,:,kc,ic) * t_di(:,:,:,:,kc,jc) - bxli(:,:,:,:,kc,ic) * t_dr(:,:,:,:,kc,jc)

               end do
            end do
         end do

         sxlr = 0.0d0
         sxli = 0.0d0
         sbxlr = cshift(bxlr,dim=that,shift=it)
         sbxli = cshift(bxli,dim=that,shift=it)

         do ic = 1,nc
            do jc = 1,nc
               do kc = 1,nc

                  sxlr(:,:,:,:,ic,jc) = sxlr(:,:,:,:,ic,jc)
     &                 + tmp1r(:,:,:,:,ic,kc) * sbxlr(:,:,:,:,kc,jc) - tmp1i(:,:,:,:,ic,kc) * sbxli(:,:,:,:,kc,jc)

                  sxli(:,:,:,:,ic,jc) = sxli(:,:,:,:,ic,jc)
     &                 + tmp1r(:,:,:,:,ic,kc) * sbxli(:,:,:,:,kc,jc) + tmp1i(:,:,:,:,ic,kc) * sbxlr(:,:,:,:,kc,jc)

               end do
            end do
         end do
