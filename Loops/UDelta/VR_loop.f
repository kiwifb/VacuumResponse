!$omp parallel setcions
!$omp section
      call delta_loops(xhat,yhat,utr,uti,WTp_xy(:,:,:,:,:,:))
!$omp section
      call delta_loops(xhat,zhat,utr,uti,WTp_xz(:,:,:,:,:,:))
!$omp end parallel setcions

!$omp parallel do
      do ilp = 0, nloop

         do it = 1, nL(that)

            WTpavg_xy(ilp,it) = Sum( WTp_xy(:,:,:,:,ilp,it) ) / volume
            WTpavg_xz(ilp,it) = Sum( WTp_xz(:,:,:,:,ilp,it) ) / volume

         end do              ! it = 2, nL(that)  loop end
      end do
!$omp end parallel do
