      do ilp = 1, nloop

         actionC(iy,iz,it,ilp,Tee,offset) = sum( WY(:,:,:,:,ilp,Tee) * actionTz(:,:,:,:) ) / volume

         topChgC(iy,iz,it,ilp,Tee,offset) = sum( WY(:,:,:,:,ilp,Tee) * topChgTz(:,:,:,:) ) / volume

      end do     ! ilp loop end
