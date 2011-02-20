      do ilp = 1, nloop

         actionTC(iy,iz,it,ilp,Tee,offset) = sum( WT(:,:,:,:,ilp,Tee) * actionTz(:,:,:,:) ) / volume

         topChgTC(iy,iz,it,ilp,Tee,offset) = sum( WT(:,:,:,:,ilp,Tee) * topChgTz(:,:,:,:) ) / volume

      end do     ! ilp loop end
