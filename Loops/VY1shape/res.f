         include 'leftleg.f'
         include 'xleg.f'

         Res = 0.0d0
         do ic = 1,36

            Res(:,:,:,:) = Res(:,:,:,:)  + fper(ic) * (
     &             t_ur(:,:,:,:,la(ic),lap(ic)) * sllr(:,:,:,:,lb(ic),lbp(ic)) * sxlr(:,:,:,:,lc(ic),lcp(ic))
     &           - t_ur(:,:,:,:,la(ic),lap(ic)) * slli(:,:,:,:,lb(ic),lbp(ic)) * sxli(:,:,:,:,lc(ic),lcp(ic))
     &           - t_ui(:,:,:,:,la(ic),lap(ic)) * sllr(:,:,:,:,lb(ic),lbp(ic)) * sxli(:,:,:,:,lc(ic),lcp(ic))
     &           - t_ui(:,:,:,:,la(ic),lap(ic)) * slli(:,:,:,:,lb(ic),lbp(ic)) * sxlr(:,:,:,:,lc(ic),lcp(ic))
     &           )

         end do
