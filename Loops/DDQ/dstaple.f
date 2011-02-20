         tmpr = 0.0d0
         tmpi = 0.0d0

         do ic = 1,nc
            do jc = 1,nc
               do kc = 1,nc

                  tmpr(:,:,:,:,ic,jc) = tmpr(:,:,:,:,ic,jc)  &
                       + llr(:,:,:,:,kc,ic) * stlr(:,:,:,:,kc,jc) + lli(:,:,:,:,kc,ic) * stli(:,:,:,:,kc,jc)

                  tmpi(:,:,:,:,ic,jc) = tmpi(:,:,:,:,ic,jc)   &
                       + llr(:,:,:,:,kc,ic) * stli(:,:,:,:,kc,jc) - lli(:,:,:,:,kc,ic) * stlr(:,:,:,:,kc,jc)

               end do
            end do
         end do

         str = 0.0d0
         sti = 0.0d0
         slr = cshift(llr, dim= that, shift= it)
         sli = cshift(lli, dim= that, shift= it)

         do ic = 1,nc
            do jc = 1,nc
               do kc = 1,nc

                  str(:,:,:,:,ic,jc) = str(:,:,:,:,ic,jc)  &
                       + tmpr(:,:,:,:,ic,kc) * slr(:,:,:,:,kc,jc) - tmpi(:,:,:,:,ic,kc) * sli(:,:,:,:,kc,jc)

                  sti(:,:,:,:,ic,jc) = sti(:,:,:,:,ic,jc)  &
                       + tmpr(:,:,:,:,ic,kc) * sli(:,:,:,:,kc,jc) + tmpi(:,:,:,:,ic,kc) * slr(:,:,:,:,kc,jc)

               end do
            end do
         end do
