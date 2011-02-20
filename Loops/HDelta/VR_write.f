!$omp parallel do
      do iy = 0, nz-1
         do iz = 0, ny-1

            icy = modulo(nz - iz , nz)
            icz = modulo(     iy , ny)

            actionTpC_xy(iz,iy,:,1:nloop,:,:) = actionTpC_xy(iz,iy,:,1:nloop,:,:) + actionTpC_xz(icz,icy,:,1:nloop,:,:)

            topChgTpC_xy(iz,iy,:,1:nloop,:,:) = topChgTpC_xy(iz,iy,:,1:nloop,:,:) + topChgTpC_xz(icz,icy,:,1:nloop,:,:)

         enddo
      enddo
!$omp end parallel do
c
!$omp parallel workshare
      WTpavg_xy(1:nloop,:) = WTpavg_xy(1:nloop,:) + WTpavg_xz(1:nloop,:)
!$omp end parallel workshare
      filename=trim(output)//'.DLTavg4'//trim(fstr1)
      call writeshape(filename,actionTpC_xy(:,:,:,:,:,:),WTpavg_xy(:,:),avgAction,offmax)

      filename=trim(output)//'.DLTavg4'//trim(fstr2)
      call writeshape(filename,TopChgTpC_xy(:,:,:,:,:,:),WTpavg_xy(:,:),avgTopChg,offmax)
