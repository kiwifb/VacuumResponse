c
c        kink 1
c
         st1r = 0.0d0
         st1i = 0.0d0

         do ic = 1,nc
            do jc = 1,nc
               do kc = 1,nc

                  st1r(:,:,:,:,ic,jc) = st1r(:,:,:,:,ic,jc)
     &                 + t_1r(:,:,:,:,ic,kc) * sbllr(:,:,:,:,jc,kc) + t_1i(:,:,:,:,ic,kc) * sblli(:,:,:,:,jc,kc)

                  st1i(:,:,:,:,ic,jc) = st1i(:,:,:,:,ic,jc)
     &                 - t_1r(:,:,:,:,ic,kc) * sblli(:,:,:,:,jc,kc) + t_1i(:,:,:,:,ic,kc) * sbllr(:,:,:,:,jc,kc)

               end do
            end do
         end do
c
c        kink 2
c
         st2r = 0.0d0
         st2i = 0.0d0

         do ic = 1,nc
            do jc = 1,nc
               do kc = 1,nc

                  st2r(:,:,:,:,ic,jc) = st2r(:,:,:,:,ic,jc)
     &                 + bllr(:,:,:,:,kc,ic) * t_2r(:,:,:,:,kc,jc) + blli(:,:,:,:,kc,ic) * t_2i(:,:,:,:,kc,jc)

                  st2i(:,:,:,:,ic,jc) = st2i(:,:,:,:,ic,jc)
     &                 - blli(:,:,:,:,kc,ic) * t_2r(:,:,:,:,kc,jc) + bllr(:,:,:,:,kc,ic) * t_2i(:,:,:,:,kc,jc)

               end do
            end do
         end do
c
c        kinked staple 3
c
         tmp1r = 0.0d0
         tmp1i = 0.0d0

         do ic = 1,nc
            do jc = 1,nc
               do kc = 1,nc

                  tmp1r(:,:,:,:,ic,jc) = tmp1r(:,:,:,:,ic,jc)
     &                 + bxlr(:,:,:,:,kc,ic) * t_3r(:,:,:,:,kc,jc) + bxli(:,:,:,:,kc,ic) * t_3i(:,:,:,:,kc,jc)

                  tmp1i(:,:,:,:,ic,jc) = tmp1i(:,:,:,:,ic,jc)
     &                 + bxlr(:,:,:,:,kc,ic) * t_3i(:,:,:,:,kc,jc) - bxli(:,:,:,:,kc,ic) * t_3r(:,:,:,:,kc,jc)

               end do
            end do
         end do

         st3r = 0.0d0
         st3i = 0.0d0

         do ic = 1,nc
            do jc = 1,nc
               do kc = 1,nc

                  st3r(:,:,:,:,ic,jc) = st3r(:,:,:,:,ic,jc)
     &                 + tmp1r(:,:,:,:,ic,kc) * sbrlr(:,:,:,:,jc,kc) + tmp1i(:,:,:,:,ic,kc) * sbrli(:,:,:,:,jc,kc)

                  st3i(:,:,:,:,ic,jc) = st3i(:,:,:,:,ic,jc)
     &                 - tmp1r(:,:,:,:,ic,kc) * sbrli(:,:,:,:,jc,kc) + tmp1i(:,:,:,:,ic,kc) * sbrlr(:,:,:,:,jc,kc)

               end do
            end do
         end do
c
c        loop
c
         do ic = 1,36

            Res(:,:,:,:) = Res(:,:,:,:)  + fper(ic) * (
     &             st2r(:,:,:,:,la(ic),lap(ic)) * st1r(:,:,:,:,lb(ic),lbp(ic)) * st3r(:,:,:,:,lc(ic),lcp(ic))
     &           - st2r(:,:,:,:,la(ic),lap(ic)) * st1i(:,:,:,:,lb(ic),lbp(ic)) * st3i(:,:,:,:,lc(ic),lcp(ic))
     &           - st2i(:,:,:,:,la(ic),lap(ic)) * st1r(:,:,:,:,lb(ic),lbp(ic)) * st3i(:,:,:,:,lc(ic),lcp(ic))
     &           - st2i(:,:,:,:,la(ic),lap(ic)) * st1i(:,:,:,:,lb(ic),lbp(ic)) * st3r(:,:,:,:,lc(ic),lcp(ic))
     &           )

         end do

