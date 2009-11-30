
      actionLC(:,:,:,1:nloop,:,:) = actionLC(:,:,:,1:nloop,:,:) / (8 * volume)
!      topChgLC(:,:,:,1:nloop,:,:) = topChgLC(:,:,:,1:nloop,:,:) / (8 * volume)

      write(filename,'(a,a,a)') trim(output),'.LPavg8',trim(fstr1)
      call writeshape(filename,actionLC(:,:,:,:,:,:),WLavg(:,:),avgAction,offmax)

!      write(filename,'(a,a,a)') trim(output),'.LPavg8',trim(fstr2)
!      call writeshape(filename,TopChgLC(:,:,:,:,:,:),WLavg(:,:),avgTopChg,offmax)
