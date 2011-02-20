      do ilp = 1, nloop

        actionTOP_C(iy,iz,it,ilp,Tee,offset) = sum(WOP(:,:,:,:,ilp,Tee) * actionTz(:,:,:,:) ) / volume

c        topChgTOP_C(iy,iz,it,ilp,Tee,offset) = sum(WOP(:,:,:,:,ilp,Tee) * topChgTz(:,:,:,:) ) / volume

      end do ! idq
