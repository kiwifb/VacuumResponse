c
c  Time averaged Quantities
c
      double precision,dimension(nx,ny,nz,nt,1:ysep,1:dqsep,nL(that)) :: WPLxy,WPLxz ! planar wilson loops
      double precision,dimension(nx,ny,nz,nt,1:ysep,1:dqsep,nL(that)) :: WOP ! non-planar wilson loops
!HPF$ ALIGN (:,:,*,*,*,*,*) WITH wloops(:,:,*,*,*) :: WPLxy,WPLxz,WOP
      double precision,dimension(1:ysep,1:dqsep,nL(that))             :: WPLxy_avg,WPLxz_avg,WOP_avg
!HPF$ DISTRIBUTE WPLxy_avg(BLOCK,BLOCK,*)
!HPF$ DISTRIBUTE WPLxz_avg(BLOCK,BLOCK,*)
!HPF$ DISTRIBUTE WOP_avg(BLOCK,BLOCK,*)
c
c  Correlations
c
      double precision,dimension(0:ny-1,0:nz-1,0:nt-1,1:ysep,1:dqsep,nL(that),offmax) :: actionTPLxy_C,actionTPLxz_C
      double precision,dimension(0:ny-1,0:nz-1,0:nt-1,1:ysep,1:dqsep,nL(that),offmax) :: actionTOP_C
      double precision,dimension(0:ny-1,0:nz-1,0:nt-1,1:ysep,1:dqsep,nL(that),offmax) :: topChgTPLxy_C,topChgTPLxz_C
      double precision,dimension(0:ny-1,0:nz-1,0:nt-1,1:ysep,1:dqsep,nL(that),offmax) :: topChgTOP_C
!HPF$ ALIGN (*,:,:,*,*,*,*) WITH correl(*,:,:,*,*) :: actionTPLxy_C,topChgTPLxy_C,actionTPLxz_C,topChgTPLxz_C
!HPF$ ALIGN (*,:,:,*,*,*,*) WITH correl(*,:,:,*,*) :: actionTOP_C,topChgTOP_C
c
c  Counters
c
      integer                                                           :: idq, iqq

