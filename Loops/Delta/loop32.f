c
c        kink 1
c
         st1r = 0.0d0
         st1i = 0.0d0

         do ic = 1,nc
            do jc = 1,nc
               do kc = 1,nc

                  st1r(:,:,:,:,ic,jc) = st1r(:,:,:,:,ic,jc)
     &                 + bxlr(:,:,:,:,ic,kc) * t_1r(:,:,:,:,kc,jc) - bxli(:,:,:,:,ic,kc) * t_1i(:,:,:,:,kc,jc)

                  st1i(:,:,:,:,ic,jc) = st1i(:,:,:,:,ic,jc)
     &                 + bxli(:,:,:,:,ic,kc) * t_1r(:,:,:,:,kc,jc) + bxlr(:,:,:,:,ic,kc) * t_1i(:,:,:,:,kc,jc)

               end do
            end do
         end do
c
c        kinked staple 2
c
         tmp1r = 0.0d0
         tmp1i = 0.0d0

         do ic = 1,nc
            do jc = 1,nc
               do kc = 1,nc

                  tmp1r(:,:,:,:,ic,jc) = tmp1r(:,:,:,:,ic,jc)
     &                 + brlr(:,:,:,:,kc,ic) * t_2r(:,:,:,:,kc,jc) + brli(:,:,:,:,kc,ic) * t_2i(:,:,:,:,kc,jc)

                  tmp1i(:,:,:,:,ic,jc) = tmp1i(:,:,:,:,ic,jc)
     &                 - brli(:,:,:,:,kc,ic) * t_2r(:,:,:,:,kc,jc) + brlr(:,:,:,:,kc,ic) * t_2i(:,:,:,:,kc,jc)

               end do
            end do
         end do

         st2r = 0.0d0
         st2i = 0.0d0

         do ic = 1,nc
            do jc = 1,nc
               do kc = 1,nc

                  st2r(:,:,:,:,ic,jc) = st2r(:,:,:,:,ic,jc)
     &                 + tmp1r(:,:,:,:,ic,kc) * sbllr(:,:,:,:,kc,jc) - tmp1i(:,:,:,:,ic,kc) * sblli(:,:,:,:,kc,jc)

                  st2i(:,:,:,:,ic,jc) = st2i(:,:,:,:,ic,jc)
     &                 + tmp1r(:,:,:,:,ic,kc) * sblli(:,:,:,:,kc,jc) + tmp1i(:,:,:,:,ic,kc) * sbllr(:,:,:,:,kc,jc)

               end do
            end do
         end do
c
c        kink 3
c
         st3r = 0.0d0
         st3i = 0.0d0

         do ic = 1,nc
            do jc = 1,nc
               do kc = 1,nc

                  st3r(:,:,:,:,ic,jc) = st3r(:,:,:,:,ic,jc)
     &                 + t_3r(:,:,:,:,ic,kc) * sbxlr(:,:,:,:,kc,jc) - t_3i(:,:,:,:,ic,kc) * sbxli(:,:,:,:,kc,jc)

                  st3i(:,:,:,:,ic,jc) = st3i(:,:,:,:,ic,jc)
     &                 + t_3i(:,:,:,:,ic,kc) * sbxlr(:,:,:,:,kc,jc) + t_3r(:,:,:,:,ic,kc) * sbxli(:,:,:,:,kc,jc)

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
