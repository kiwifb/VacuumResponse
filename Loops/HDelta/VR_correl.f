      do ilp = 1, nloop

         actionTpC_xy(iy,iz,it,ilp,Tee,offset) = sum( WTp_xy(:,:,:,:,ilp,Tee) * actionTz(:,:,:,:) ) / volume
         actionTpC_xz(iy,iz,it,ilp,Tee,offset) = sum( WTp_xz(:,:,:,:,ilp,Tee) * actionTz(:,:,:,:) ) / volume

         topChgTpC_xy(iy,iz,it,ilp,Tee,offset) = sum( WTp_xy(:,:,:,:,ilp,Tee) * topChgTz(:,:,:,:) ) / volume
         topChgTpC_xz(iy,iz,it,ilp,Tee,offset) = sum( WTp_xz(:,:,:,:,ilp,Tee) * topChgTz(:,:,:,:) ) / volume

      end do     ! iylp loop end
