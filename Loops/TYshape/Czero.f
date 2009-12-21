c
c     ------------------------------
c     CASE 0 - three quarks at (0,0)
c     ------------------------------
c
         Res = 0.0d0
         do ic = 1,36

            Res(:,:,:,:) = Res(:,:,:,:)  + fper(ic) * (
     &             uptr(:,:,:,:,la(ic),lap(ic)) * uptr(:,:,:,:,lb(ic),lbp(ic)) * uptr(:,:,:,:,lc(ic),lcp(ic))
     &           - uptr(:,:,:,:,la(ic),lap(ic)) * upti(:,:,:,:,lb(ic),lbp(ic)) * upti(:,:,:,:,lc(ic),lcp(ic))
     &           - upti(:,:,:,:,la(ic),lap(ic)) * uptr(:,:,:,:,lb(ic),lbp(ic)) * upti(:,:,:,:,lc(ic),lcp(ic))
     &           - upti(:,:,:,:,la(ic),lap(ic)) * upti(:,:,:,:,lb(ic),lbp(ic)) * uptr(:,:,:,:,lc(ic),lcp(ic))
     &           )

         end do

