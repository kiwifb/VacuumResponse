!$omp parallel do
      do ix = 0, nt-1

         icx = modulo(nt - ix , nt)

         actionTp_xy_C(:,:,ix,1:nloop,:,:) = actionTp_xy_C(:,:, ix,1:nloop,:,:) + actionTm_xy_C(:,:,icx,1:nloop,:,:)

         topChgTp_xy_C(:,:,ix,1:nloop,:,:) = topChgTp_xy_C(:,:, ix,1:nloop,:,:) + topChgTm_xy_C(:,:,icx,1:nloop,:,:)

         actionT4p_xy_C(:,:,ix,1:nloop,:,:) = actionT4p_xy_C(:,:, ix,1:nloop,:,:) + actionT4m_xy_C(:,:,icx,1:nloop,:,:)

         topChgT4p_xy_C(:,:,ix,1:nloop,:,:) = topChgT4p_xy_C(:,:, ix,1:nloop,:,:) + topChgT4m_xy_C(:,:,icx,1:nloop,:,:)

      end do
!$omp end parallel do
c
!$omp parallel setcions
!$omp section
      WTp_xy_avg(1:nloop,:) = WTp_xy_avg(1:nloop,:) + WTm_xy_avg(1:nloop,:)
!$omp section
      WT4p_xy_avg(1:nloop,:) = WT4p_xy_avg(1:nloop,:) + WT4m_xy_avg(1:nloop,:)
!$omp end parallel setcions

!$omp parallel do
      do iy = 0, nz-1
         do iz = 0, ny-1

            icy = modulo(nz - iz , nz)
            icz = modulo(     iy , ny)

            actionTp_xy_C(iz,iy,:,1:nloop,:,:) = actionTp_xy_C(iz,iy,:,1:nloop,:,:) + actionTp_xz_C(icz,icy,:,1:nloop,:,:)

            topChgTp_xy_C(iz,iy,:,1:nloop,:,:) = topChgTp_xy_C(iz,iy,:,1:nloop,:,:) + topChgTp_xz_C(icz,icy,:,1:nloop,:,:)

            actionT4p_xy_C(iz,iy,:,1:nloop,:,:) = actionT4p_xy_C(iz,iy,:,1:nloop,:,:) + actionT4p_xz_C(icz,icy,:,1:nloop,:,:)

            topChgT4p_xy_C(iz,iy,:,1:nloop,:,:) = topChgT4p_xy_C(iz,iy,:,1:nloop,:,:) + topChgT4p_xz_C(icz,icy,:,1:nloop,:,:)

         end do
      end do
!$omp end parallel do
c
!$omp parallel setcions
!$omp section
      WTp_xy_avg(1:nloop,:) = WTp_xy_avg(1:nloop,:) + WTp_xz_avg(1:nloop,:)
!$omp section
      WT4p_xy_avg(1:nloop,:) = WT4p_xy_avg(1:nloop,:) + WT4p_xz_avg(1:nloop,:)
!$omp end parallel setcions

!$omp parallel do
      do iy = 0, nz-1
         do iz = 0, ny-1
            do ix = 0, nt-1

               icx = modulo(nt - ix , nt)
               icy = modulo(nz - iz , nz)
               icz = modulo(     iy , ny)

               actionTp_xy_C(iz,iy,ix,1:nloop,:,:) = actionTp_xy_C(iz,iy,ix,1:nloop,:,:) + actionTm_xz_C(icz,icy,icx,1:nloop,:,:)

               topChgTp_xy_C(iz,iy,ix,1:nloop,:,:) = topChgTp_xy_C(iz,iy,ix,1:nloop,:,:) + topChgTm_xz_C(icz,icy,icx,1:nloop,:,:)

               actionT4p_xy_C(iz,iy,ix,1:nloop,:,:) = actionT4p_xy_C(iz,iy,ix,1:nloop,:,:) + actionT4m_xz_C(icz,icy,icx,1:nloop,:,:)

               topChgT4p_xy_C(iz,iy,ix,1:nloop,:,:) = topChgT4p_xy_C(iz,iy,ix,1:nloop,:,:) + topChgT4m_xz_C(icz,icy,icx,1:nloop,:,:)

            end do
         end do
      end do
!$omp end parallel do
c
!$omp parallel workshare
      WTp_xy_avg(1:nloop,:) = WTp_xy_avg(1:nloop,:) + WTm_xz_avg(1:nloop,:)
      WT4p_xy_avg(1:nloop,:) = WT4p_xy_avg(1:nloop,:) + WT4m_xz_avg(1:nloop,:)

      actionTp_xy_C(:,:,:,1:nloop,:,:) = actionTp_xy_C(:,:,:,1:nloop,:,:) / 4
      topChgTp_xy_C(:,:,:,1:nloop,:,:) = topChgTp_xy_C(:,:,:,1:nloop,:,:) / 4
      WTp_xy_avg(1:nloop,:) = WTp_xy_avg(1:nloop,:) / 4

      actionT4p_xy_C(:,:,:,1:nloop,:,:) = actionT4p_xy_C(:,:,:,1:nloop,:,:) / 4
      topChgT4p_xy_C(:,:,:,1:nloop,:,:) = topChgT4p_xy_C(:,:,:,1:nloop,:,:) / 4
      WT4p_xy_avg(1:nloop,:) = WT4p_xy_avg(1:nloop,:) / 4

      actionTp_xy_C = cshift(cshift(actionTp_xy_C,dim=3,shift= 4),dim=2,shift= ny/2)
      actionT4p_xy_C = cshift(cshift(actionT4p_xy_C,dim=3,shift= 4),dim=2,shift= ny/2)

      topChgTp_xy_C = cshift(cshift(topChgTp_xy_C,dim=3,shift= 4),dim=2,shift= ny/2)
      topChgT4p_xy_C = cshift(cshift(topChgT4p_xy_C,dim=3,shift= 4),dim=2,shift= ny/2)
!$omp end parallel workshare
      filename=trim(output)//'.DQ2avg4'//trim(fstr1)
      call writeshape(filename,actionTp_xy_C(:,:,:,:,:,:),WTp_xy_avg(:,:),avgAction,offmax)

      filename=trim(output)//'.DQ2avg4'//trim(fstr2)
      call writeshape(filename,TopChgTp_xy_C(:,:,:,:,:,:),WTp_xy_avg(:,:),avgTopChg,offmax)

      filename=trim(output)//'.DQ4avg4'//trim(fstr1)
      call writeshape(filename,actionT4p_xy_C(:,:,:,:,:,:),WT4p_xy_avg(:,:),avgAction,offmax)

      filename=trim(output)//'.DQ4avg4'//trim(fstr2)
      call writeshape(filename,TopChgT4p_xy_C(:,:,:,:,:,:),WT4p_xy_avg(:,:),avgTopChg,offmax)
