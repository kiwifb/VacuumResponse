!$omp parallel setcions
!$omp section
      call t_loops(xhat,yhat,utr,uti,
     &             WTp_xy(:,:,:,:,:,:),WTm_xy(:,:,:,:,:,:),WT4p_xy(:,:,:,:,:,:),WT4m_xy(:,:,:,:,:,:))
!$omp section
      call t_loops(xhat,zhat,utr,uti,
     &             WTp_xz(:,:,:,:,:,:),WTm_xz(:,:,:,:,:,:),WT4p_xz(:,:,:,:,:,:),WT4m_xz(:,:,:,:,:,:))
!$omp end parallel setcions

!$omp parallel do
      do ilp = 1, nloop

         do it = 1, nL(that)

            WTp_xy_avg(ilp,it)  = Sum( WTp_xy(:,:,:,:,ilp,it) )  / volume
            WTm_xy_avg(ilp,it)  = Sum( WTm_xy(:,:,:,:,ilp,it) )  / volume
            WTp_xz_avg(ilp,it)  = Sum( WTp_xz(:,:,:,:,ilp,it) )  / volume
            WTm_xz_avg(ilp,it)  = Sum( WTm_xz(:,:,:,:,ilp,it) )  / volume
            WT4p_xy_avg(ilp,it) = Sum( WT4p_xy(:,:,:,:,ilp,it) ) / volume
            WT4m_xy_avg(ilp,it) = Sum( WT4m_xy(:,:,:,:,ilp,it) ) / volume
            WT4p_xz_avg(ilp,it) = Sum( WT4p_xz(:,:,:,:,ilp,it) ) / volume
            WT4m_xz_avg(ilp,it) = Sum( WT4m_xz(:,:,:,:,ilp,it) ) / volume

         end do              ! it = 2, nL(that)  loop end
      end do
!$omp end parallel do
