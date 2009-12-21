c
c        staple 1
c
         tmp1r = 0.0d0
         tmp1i = 0.0d0

         do ic = 1,nc
            do jc = 1,nc
               do kc = 1,nc

                  tmp1r(:,:,:,:,ic,jc) = tmp1r(:,:,:,:,ic,jc)
     &                 + bxlr(:,:,:,:,ic,kc) * t_1r(:,:,:,:,kc,jc) - bxli(:,:,:,:,ic,kc) * t_1i(:,:,:,:,kc,jc)

                  tmp1i(:,:,:,:,ic,jc) = tmp1i(:,:,:,:,ic,jc)
     &                 + bxlr(:,:,:,:,ic,kc) * t_1i(:,:,:,:,kc,jc) + bxli(:,:,:,:,ic,kc) * t_1r(:,:,:,:,kc,jc)

               end do
            end do
         end do

         st1r = 0.0d0
         st1i = 0.0d0

         do ic = 1,nc
            do jc = 1,nc
               do kc = 1,nc

                  st1r(:,:,:,:,ic,jc) = st1r(:,:,:,:,ic,jc)
     &                 + tmp1r(:,:,:,:,ic,kc) * sbxlr(:,:,:,:,jc,kc) + tmp1i(:,:,:,:,ic,kc) * sbxli(:,:,:,:,jc,kc)

                  st1i(:,:,:,:,ic,jc) = st1i(:,:,:,:,ic,jc)
     &                 - tmp1r(:,:,:,:,ic,kc) * sbxli(:,:,:,:,jc,kc) + tmp1i(:,:,:,:,ic,kc) * sbxlr(:,:,:,:,jc,kc)

               end do
            end do
         end do
c
c        staple 2
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
     &                 + tmp1r(:,:,:,:,ic,kc) * sbrlr(:,:,:,:,kc,jc) - tmp1i(:,:,:,:,ic,kc) * sbrli(:,:,:,:,kc,jc)

                  st2i(:,:,:,:,ic,jc) = st2i(:,:,:,:,ic,jc)
     &                 + tmp1r(:,:,:,:,ic,kc) * sbrli(:,:,:,:,kc,jc) + tmp1i(:,:,:,:,ic,kc) * sbrlr(:,:,:,:,kc,jc)

               end do
            end do
         end do
c
c        loop
c
         do ic = 1,36

            Res(:,:,:,:) = Res(:,:,:,:)  + fper(ic) * (
     &             t_3r(:,:,:,:,la(ic),lap(ic)) * st1r(:,:,:,:,lb(ic),lbp(ic)) * st2r(:,:,:,:,lc(ic),lcp(ic))
     &           - t_3r(:,:,:,:,la(ic),lap(ic)) * st1i(:,:,:,:,lb(ic),lbp(ic)) * st2i(:,:,:,:,lc(ic),lcp(ic))
     &           - t_3i(:,:,:,:,la(ic),lap(ic)) * st1r(:,:,:,:,lb(ic),lbp(ic)) * st2i(:,:,:,:,lc(ic),lcp(ic))
     &           - t_3i(:,:,:,:,la(ic),lap(ic)) * st1i(:,:,:,:,lb(ic),lbp(ic)) * st2r(:,:,:,:,lc(ic),lcp(ic))
     &           )

         end do
