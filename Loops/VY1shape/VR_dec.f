      double precision,dimension(nx,ny,nz,nt,0:nloop,nL(that))   :: WY ! wilson loops
!HPF$ ALIGN (:,:,*,*,*,*) WITH wloops(:,:,*,*,*) :: WY
      double precision,dimension(0:nloop,nL(that))               :: WYavg
!HPF$ DISTRIBUTE WYavg(BLOCK,*)
c
c  Correlations
c
      double precision,dimension(0:ny-1,0:nz-1,0:nt-1,0:nloop,nL(that),offmax) :: actionC
      double precision,dimension(0:ny-1,0:nz-1,0:nt-1,0:nloop,nL(that),offmax) :: topChgC
!HPF$ ALIGN (*,:,:,*,*,*) WITH correl(*,:,:,*,*) :: actionC,topChgC
