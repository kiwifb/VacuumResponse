      call quad_loops(xhat,yhat,utr,uti,WT)

!$omp parallel do
      do ilp = 1, nloop

         do it = 1, nL(that)

            WTavg(ilp,it)  = Sum( WT(:,:,:,:,ilp,it) )  / volume

         end do              ! it = 2, nL(that)  loop end
      end do
!$omp end parallel do
