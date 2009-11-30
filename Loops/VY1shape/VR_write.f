      filename=trim(output)//'.VY1avg4'//trim(fstr1)
      call writeshape(filename,actionC(:,:,:,:,:,:),WYavg(:,:),avgAction,offmax)

      filename=trim(output)//'.VY1avg4'//trim(fstr2)
      call writeshape(filename,TopChgC(:,:,:,:,:,:),WYavg(:,:),avgTopChg,offmax)
