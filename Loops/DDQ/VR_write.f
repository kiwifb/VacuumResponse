
      filename=trim(output)//'.DDQ-PL'//trim(fstr1)
      call writeshape(filename,actionTPLxy_C(:,:,:,:,:,:,:),WPLxy_avg(:,:,:),avgAction,offmax)

c      filename=trim(output)//'.DDQ-PL'//trim(fstr2)
c      call writeshape(filename,TopChgTPLxy_C(:,:,:,:,:,:,:),WPLxy_avg(:,:,:),avgTopChg,offmax)

c      filename=trim(output)//'.DDQ-OP'//trim(fstr1)
c      call writeshape(filename,actionTOP_C(:,:,:,:,:,:,:),WOP_avg(:,:,:),avgAction,offmax)

c      filename=trim(output)//'.DQ4avg4'//trim(fstr2)
c      call writeshape(filename,TopChgTOP_C(:,:,:,:,:,:,:),WOP_avg(:,:,:),avgTopChg,offmax)
