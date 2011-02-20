      double precision,dimension(nx,ny,nz,nt,0:nloop,nL(that))   :: WT ! wilson loops
!HPF$ ALIGN (:,:,*,*,*,*) WITH wloops(:,:,*,*,*) :: WT
      double precision,dimension(0:nloop,nL(that))              :: WTavg
!HPF$ DISTRIBUTE WTavg(BLOCK,*)
c
c  Correlations
c
      double precision,dimension(0:ny-1,0:nz-1,0:nt-1,0:nloop,nL(that),offmax) :: actionTC
      double precision,dimension(0:ny-1,0:nz-1,0:nt-1,0:nloop,nL(that),offmax) :: topChgTC
!HPF$ ALIGN (*,:,:,*,*,*) WITH correl(*,:,:,*,*) :: actionTC,topChgTC

