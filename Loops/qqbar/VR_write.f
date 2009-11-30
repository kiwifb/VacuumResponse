!$omp parallel workshare
      actionTC = cshift(cshift(actionTC,dim=3,shift= 4),dim=2,shift= ny/2)

      topChgTC = cshift(cshift(topChgTC,dim=3,shift= 4),dim=2,shift= ny/2)
!$omp end parallel workshare

      filename=trim(output)//'.QQ'//trim(fstr1)
      call writeshape(filename,actionTC(:,:,:,:,:,:),WTavg(:,:),avgAction,offmax)

      filename=trim(output)//'.QQ'//trim(fstr2)
      call writeshape(filename,TopChgTC(:,:,:,:,:,:),WTavg(:,:),avgTopChg,offmax)

