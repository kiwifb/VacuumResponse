      call vy_loops(xhat,yhat,ur,ui,WY(:,:,:,:,:,:))

      do ilp = 0, nloop
         do it = 1, nL(that)

            WYavg(ilp,it) = Sum( WY(:,:,:,:,ilp,it) ) / volume

         end do              ! it = 2, nL(that)  loop end
      end do
