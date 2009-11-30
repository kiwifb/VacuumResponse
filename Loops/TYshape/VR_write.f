!$omp parallel do
      do ix = 0, nt-1

         icx = modulo(nt - ix , nt)

         actionYp_xy_C(:,:,ix,1:nloop,:,:) = actionYp_xy_C(:,:, ix,1:nloop,:,:) + actionYm_xy_C(:,:,icx,1:nloop,:,:)

         topChgYp_xy_C(:,:,ix,1:nloop,:,:) = topChgYp_xy_C(:,:, ix,1:nloop,:,:) + topChgYm_xy_C(:,:,icx,1:nloop,:,:)

      enddo
!$omp end parallel do
c
!$omp parallel workshare
       WYp_xy_avg(1:nloop,:) = WYp_xy_avg(1:nloop,:) + WYm_xy_avg(1:nloop,:)
!$omp end parallel workshare
 
!$omp parallel do
      do iy = 0, nz-1
         do iz = 0, ny-1

            icy = modulo(nz - iz , nz)
            icz = modulo(     iy , ny)

            actionYp_xy_C(iz,iy,:,1:nloop,:,:) = actionYp_xy_C(iz,iy,:,1:nloop,:,:) + actionYp_xz_C(icz,icy,:,1:nloop,:,:)

            topChgYp_xy_C(iz,iy,:,1:nloop,:,:) = topChgYp_xy_C(iz,iy,:,1:nloop,:,:) + topChgYp_xz_C(icz,icy,:,1:nloop,:,:)

         enddo
      enddo
!$omp parallel do
c
!$omp parallel workshare
       WYp_xy_avg(1:nloop,:) = WYp_xy_avg(1:nloop,:) + WYp_xz_avg(1:nloop,:)
!$omp end parallel workshare
 
!$omp parallel do
      do iy = 0, nz-1
         do iz = 0, ny-1
            do ix = 0, nt-1

               icx = modulo(nt - ix , nt)
               icy = modulo(nz - iz , nz)
               icz = modulo(     iy , ny)

               actionYp_xy_C(iz,iy,ix,1:nloop,:,:) = actionYp_xy_C(iz,iy,ix,1:nloop,:,:) + actionYm_xz_C(icz,icy,icx,1:nloop,:,:)

               topChgYp_xy_C(iz,iy,ix,1:nloop,:,:) = topChgYp_xy_C(iz,iy,ix,1:nloop,:,:) + topChgYm_xz_C(icz,icy,icx,1:nloop,:,:)

            enddo
         enddo
      enddo
!$omp parallel do
c
!$omp parallel workshare
       WYp_xy_avg(1:nloop,:) = WYp_xy_avg(1:nloop,:) + WYm_xz_avg(1:nloop,:)

      actionYp_xy_C(:,:,:,1:nloop,:,:) = actionYp_xy_C(:,:,:,1:nloop,:,:) / 4
      topChgYp_xy_C(:,:,:,1:nloop,:,:) = topChgYp_xy_C(:,:,:,1:nloop,:,:) / 4
      WYp_xy_avg(1:nloop,:) = WYp_xy_avg(1:nloop,:) / 4
!$omp end parallel workshare
 
      filename=trim(output)//'.TYavg4'//trim(fstr1)
      call writeshape(filename,actionYp_xy_C(:,:,:,:,:,:),WYp_xy_avg(:,:),avgAction,offmax)

      filename=trim(output)//'.TYavg4'//trim(fstr2)
      call writeshape(filename,TopChgYp_xy_C(:,:,:,:,:,:),WYp_xy_avg(:,:),avgTopChg,offmax)
