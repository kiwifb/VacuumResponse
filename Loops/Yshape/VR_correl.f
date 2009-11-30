      do ilp = 1, nloop

         actionYp_xy_C(iy,iz,it,ilp,Tee,offset) = sum( WYp_xy(:,:,:,:,ilp,Tee) * actionTz(:,:,:,:) ) / volume
         actionYm_xy_C(iy,iz,it,ilp,Tee,offset) = sum( WYm_xy(:,:,:,:,ilp,Tee) * actionTz(:,:,:,:) ) / volume
         actionYp_xz_C(iy,iz,it,ilp,Tee,offset) = sum( WYp_xz(:,:,:,:,ilp,Tee) * actionTz(:,:,:,:) ) / volume
         actionYm_xz_C(iy,iz,it,ilp,Tee,offset) = sum( WYm_xz(:,:,:,:,ilp,Tee) * actionTz(:,:,:,:) ) / volume

         topChgYp_xy_C(iy,iz,it,ilp,Tee,offset) = sum( WYp_xy(:,:,:,:,ilp,Tee) * topChgTz(:,:,:,:) ) / volume
         topChgYm_xy_C(iy,iz,it,ilp,Tee,offset) = sum( WYm_xy(:,:,:,:,ilp,Tee) * topChgTz(:,:,:,:) ) / volume
         topChgYp_xz_C(iy,iz,it,ilp,Tee,offset) = sum( WYp_xz(:,:,:,:,ilp,Tee) * topChgTz(:,:,:,:) ) / volume
         topChgYm_xz_C(iy,iz,it,ilp,Tee,offset) = sum( WYm_xz(:,:,:,:,ilp,Tee) * topChgTz(:,:,:,:) ) / volume

      end do     ! ilp loop end
