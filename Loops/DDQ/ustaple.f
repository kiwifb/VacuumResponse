         tmpr = 0.0d0
         tmpi = 0.0d0

         do ic = 1,nc
            do jc = 1,nc
               do kc = 1,nc

                  tmpr(:,:,:,:,ic,jc) = tmpr(:,:,:,:,ic,jc)  &
                       + llr(:,:,:,:,ic,kc) * stlr(:,:,:,:,kc,jc) - lli(:,:,:,:,ic,kc) * stli(:,:,:,:,kc,jc)

                  tmpi(:,:,:,:,ic,jc) = tmpi(:,:,:,:,ic,jc)   &
                       + llr(:,:,:,:,ic,kc) * stli(:,:,:,:,kc,jc) + lli(:,:,:,:,ic,kc) * stlr(:,:,:,:,kc,jc)

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
                       + tmpr(:,:,:,:,ic,kc) * slr(:,:,:,:,jc,kc) + tmpi(:,:,:,:,ic,kc) * sli(:,:,:,:,jc,kc)

                  sti(:,:,:,:,ic,jc) = sti(:,:,:,:,ic,jc)  &
                       - tmpr(:,:,:,:,ic,kc) * sli(:,:,:,:,jc,kc) + tmpi(:,:,:,:,ic,kc) * slr(:,:,:,:,jc,kc)

               end do
            end do
         end do
