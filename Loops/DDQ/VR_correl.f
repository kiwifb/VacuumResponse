      do idq = 1, dqsep
        do iqq = 1, ysep

          actionTPLxy_C(iy,iz,it,iqq,idq,Tee,offset) = sum( WPLxy(:,:,:,:,iqq,idq,Tee) * actionTz(:,:,:,:) ) / volume
c          actionTPLxz_C(iy,iz,it,iqq,idq,Tee,offset) = sum( WPLxz(:,:,:,:,iqq,idq,Tee) * actionTz(:,:,:,:) ) / volume
c            actionTOP_C(iy,iz,it,iqq,idq,Tee,offset) = sum(   WOP(:,:,:,:,iqq,idq,Tee) * actionTz(:,:,:,:) ) / volume

c          topChgTPLxy_C(iy,iz,it,iqq,idq,Tee,offset) = sum( WPLxy(:,:,:,:,iqq,idq,Tee) * topChgTz(:,:,:,:) ) / volume
c          topChgTPLxz_C(iy,iz,it,iqq,idq,Tee,offset) = sum( WPLxz(:,:,:,:,iqq,idq,Tee) * topChgTz(:,:,:,:) ) / volume
c            topChgTOP_C(iy,iz,it,iqq,idq,Tee,offset) = sum(   WOP(:,:,:,:,iqq,idq,Tee) * topChgTz(:,:,:,:) ) / volume

        end do ! iqq
      end do ! idq
