         include 'leftleg.f'
         include 'rightleg.f'

         Res = 0.0d0
         do ic = 1,36

            Res(:,:,:,:) = Res(:,:,:,:)  + fper(ic) * (
     &             t_tr(:,:,:,:,la(ic),lap(ic)) * sllr(:,:,:,:,lb(ic),lbp(ic)) * srlr(:,:,:,:,lc(ic),lcp(ic))
     &           - t_tr(:,:,:,:,la(ic),lap(ic)) * slli(:,:,:,:,lb(ic),lbp(ic)) * srli(:,:,:,:,lc(ic),lcp(ic))
     &           - t_ti(:,:,:,:,la(ic),lap(ic)) * sllr(:,:,:,:,lb(ic),lbp(ic)) * srli(:,:,:,:,lc(ic),lcp(ic))
     &           - t_ti(:,:,:,:,la(ic),lap(ic)) * slli(:,:,:,:,lb(ic),lbp(ic)) * srlr(:,:,:,:,lc(ic),lcp(ic))
     &           )

         end do
