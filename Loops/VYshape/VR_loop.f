!$omp parallel setcions
!$omp section
      call vy_loops_mirror(xhat,yhat,ur,ui,WYp_xy(:,:,:,:,:,:),WYm_xy(:,:,:,:,:,:))
!$omp section
      call vy_loops_mirror(xhat,zhat,ur,ui,WYp_xz(:,:,:,:,:,:),WYm_xz(:,:,:,:,:,:))
!$omp end parallel setcions

!$omp parallel do
      do ilp = 0, nloop
         do it = 1, nL(that)

            WYp_xy_avg(ilp,it) = Sum( WYp_xy(:,:,:,:,ilp,it) ) / volume
            WYm_xy_avg(ilp,it) = Sum( WYm_xy(:,:,:,:,ilp,it) ) / volume
            WYp_xz_avg(ilp,it) = Sum( WYp_xz(:,:,:,:,ilp,it) ) / volume
            WYm_xz_avg(ilp,it) = Sum( WYm_xz(:,:,:,:,ilp,it) ) / volume

         end do              ! it = 2, nL(that)  loop end
      end do
!$omp end parallel do
