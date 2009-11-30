      do ilp = 1, nloop

         actionTp_xy_C(iy,iz,it,ilp,Tee,offset) = sum( WTp_xy(:,:,:,:,ilp,Tee) * actionTz(:,:,:,:) ) / volume
         actionTm_xy_C(iy,iz,it,ilp,Tee,offset) = sum( WTm_xy(:,:,:,:,ilp,Tee) * actionTz(:,:,:,:) ) / volume
         actionTp_xz_C(iy,iz,it,ilp,Tee,offset) = sum( WTp_xz(:,:,:,:,ilp,Tee) * actionTz(:,:,:,:) ) / volume
         actionTm_xz_C(iy,iz,it,ilp,Tee,offset) = sum( WTm_xz(:,:,:,:,ilp,Tee) * actionTz(:,:,:,:) ) / volume

         topChgTp_xy_C(iy,iz,it,ilp,Tee,offset) = sum( WTp_xy(:,:,:,:,ilp,Tee) * topChgTz(:,:,:,:) ) / volume
         topChgTm_xy_C(iy,iz,it,ilp,Tee,offset) = sum( WTm_xy(:,:,:,:,ilp,Tee) * topChgTz(:,:,:,:) ) / volume
         topChgTp_xz_C(iy,iz,it,ilp,Tee,offset) = sum( WTp_xz(:,:,:,:,ilp,Tee) * topChgTz(:,:,:,:) ) / volume
         topChgTm_xz_C(iy,iz,it,ilp,Tee,offset) = sum( WTm_xz(:,:,:,:,ilp,Tee) * topChgTz(:,:,:,:) ) / volume

         actionT4p_xy_C(iy,iz,it,ilp,Tee,offset) = sum( WT4p_xy(:,:,:,:,ilp,Tee) * actionTz(:,:,:,:) ) / volume
         actionT4m_xy_C(iy,iz,it,ilp,Tee,offset) = sum( WT4m_xy(:,:,:,:,ilp,Tee) * actionTz(:,:,:,:) ) / volume
         actionT4p_xz_C(iy,iz,it,ilp,Tee,offset) = sum( WT4p_xz(:,:,:,:,ilp,Tee) * actionTz(:,:,:,:) ) / volume
         actionT4m_xz_C(iy,iz,it,ilp,Tee,offset) = sum( WT4m_xz(:,:,:,:,ilp,Tee) * actionTz(:,:,:,:) ) / volume

         topChgT4p_xy_C(iy,iz,it,ilp,Tee,offset) = sum( WT4p_xy(:,:,:,:,ilp,Tee) * topChgTz(:,:,:,:) ) / volume
         topChgT4m_xy_C(iy,iz,it,ilp,Tee,offset) = sum( WT4m_xy(:,:,:,:,ilp,Tee) * topChgTz(:,:,:,:) ) / volume
         topChgT4p_xz_C(iy,iz,it,ilp,Tee,offset) = sum( WT4p_xz(:,:,:,:,ilp,Tee) * topChgTz(:,:,:,:) ) / volume
         topChgT4m_xz_C(iy,iz,it,ilp,Tee,offset) = sum( WT4m_xz(:,:,:,:,ilp,Tee) * topChgTz(:,:,:,:) ) / volume

      end do     ! ilp loop end
