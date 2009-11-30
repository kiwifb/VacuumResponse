      do ilp = 1, nloop

         actionTpC(iy,iz,it,ilp,Tee,offset) = sum( WTp(:,:,:,:,ilp,Tee) * actionTz(:,:,:,:) ) / volume

      end do     ! iylp loop end
