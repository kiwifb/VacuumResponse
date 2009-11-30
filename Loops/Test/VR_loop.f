      call test_loops(xhat,yhat,utr,uti,WTp(:,:,:,:,:,:))

      do ilp = 0, nloop

         do it = 1, nL(that)

            WTpavg(ilp,it) = Sum( WTp(:,:,:,:,ilp,it) ) / volume

         end do              ! it = 2, nL(that)  loop end
      end do
