
      filename=trim(output)//'.4P'//trim(fstr1)
      call writeshape(filename,actionTOP_C(:,:,:,:,:,:),WOP_avg(:,:),avgAction,offmax)

c      filename=trim(output)//'.4P'//trim(fstr2)
c      call writeshape(filename,TopChgTOP_C(:,:,:,:,:,:),WOP_avg(:,:),avgTopChg,offmax)
