      double precision,dimension(nx,ny,nz,nt,0:nloop,nL(that))   :: WYp_xy,WYp_xz,WYm_xy,WYm_xz ! wilson loops
!HPF$ ALIGN (:,:,*,*,*,*) WITH wloops(:,:,*,*,*) :: WYp_xy,WYp_xz,WYm_xy,WYm_xz
      double precision,dimension(0:nloop,nL(that))               :: WYp_xy_avg,WYp_xz_avg,WYm_xy_avg,WYm_xz_avg
!HPF$ DISTRIBUTE WYp_xy_avg(BLOCK,*)
!HPF$ DISTRIBUTE WYp_xz_avg(BLOCK,*)
!HPF$ DISTRIBUTE WYm_xy_avg(BLOCK,*)
!HPF$ DISTRIBUTE WYm_xz_avg(BLOCK,*)
c
c  Correlations
c
      double precision,dimension(0:ny-1,0:nz-1,0:nt-1,0:nloop,nL(that),offmax) :: actionYp_xy_C,actionYp_xz_C
      double precision,dimension(0:ny-1,0:nz-1,0:nt-1,0:nloop,nL(that),offmax) :: actionYm_xy_C,actionYm_xz_C
      double precision,dimension(0:ny-1,0:nz-1,0:nt-1,0:nloop,nL(that),offmax) :: topChgYp_xy_C,topChgYp_xz_C
      double precision,dimension(0:ny-1,0:nz-1,0:nt-1,0:nloop,nL(that),offmax) :: topChgYm_xy_C,topChgYm_xz_C
!HPF$ ALIGN (*,:,:,*,*,*) WITH correl(*,:,:,*,*) :: actionYp_xy_C,topChgYp_xy_C
!HPF$ ALIGN (*,:,:,*,*,*) WITH correl(*,:,:,*,*) :: actionYp_xz_C,topChgYp_xz_C
!HPF$ ALIGN (*,:,:,*,*,*) WITH correl(*,:,:,*,*) :: actionYm_xy_C,topChgYm_xy_C
!HPF$ ALIGN (*,:,:,*,*,*) WITH correl(*,:,:,*,*) :: actionYm_xz_C,topChgYm_xz_C
