c
c  Time averaged Quantities
c
      double precision,dimension(nx,ny,nz,nt,0:nloop,nL(that)) :: WTp_xy,WTm_xy ! T-shape wilson loops
      double precision,dimension(nx,ny,nz,nt,0:nloop,nL(that)) :: WTp_xz,WTm_xz ! T-shape wilson loops
!HPF$ ALIGN (:,:,*,*,*,*) WITH wloops(:,:,*,*,*) :: WTp_xy,WTm_xy,WTp_xz,WTm_xz
      double precision,dimension(nx,ny,nz,nt,0:nloop,nL(that)) :: WT4p_xy,WT4m_xy ! T-shape wilson loops
      double precision,dimension(nx,ny,nz,nt,0:nloop,nL(that)) :: WT4p_xz,WT4m_xz ! T-shape wilson loops
!HPF$ ALIGN (:,:,*,*,*,*) WITH wloops(:,:,*,*,*) :: WT4p_xy,WT4m_xy,WT4p_xz,WT4m_xz

      double precision,dimension(0:nloop,nL(that))             :: WTp_xy_avg,WTm_xy_avg
      double precision,dimension(0:nloop,nL(that))             :: WTp_xz_avg,WTm_xz_avg
!HPF$ DISTRIBUTE WTp_xy_avg(BLOCK,*)
!HPF$ DISTRIBUTE WTm_xy_avg(BLOCK,*)
!HPF$ DISTRIBUTE WTp_xz_avg(BLOCK,*)
!HPF$ DISTRIBUTE WTm_xz_avg(BLOCK,*)
      double precision,dimension(0:nloop,nL(that))             :: WT4p_xy_avg,WT4m_xy_avg
      double precision,dimension(0:nloop,nL(that))             :: WT4p_xz_avg,WT4m_xz_avg
!HPF$ DISTRIBUTE WT4p_xy_avg(BLOCK,*)
!HPF$ DISTRIBUTE WT4m_xy_avg(BLOCK,*)
!HPF$ DISTRIBUTE WT4p_xz_avg(BLOCK,*)
!HPF$ DISTRIBUTE WT4m_xz_avg(BLOCK,*)
c
c  Correlations
c
      double precision,dimension(0:ny-1,0:nz-1,0:nt-1,0:nloop,nL(that),offmax) :: actionTp_xy_C,actionTm_xy_C
      double precision,dimension(0:ny-1,0:nz-1,0:nt-1,0:nloop,nL(that),offmax) :: actionTp_xz_C,actionTm_xz_C
      double precision,dimension(0:ny-1,0:nz-1,0:nt-1,0:nloop,nL(that),offmax) :: topChgTp_xy_C,topChgTm_xy_C
      double precision,dimension(0:ny-1,0:nz-1,0:nt-1,0:nloop,nL(that),offmax) :: topChgTp_xz_C,topChgTm_xz_C
!HPF$ ALIGN (*,:,:,*,*,*) WITH correl(*,:,:,*,*) :: actionTp_xy_C,topChgTp_xy_C,actionTm_xy_C,topChgTm_xy_C
!HPF$ ALIGN (*,:,:,*,*,*) WITH correl(*,:,:,*,*) :: actionTp_xz_C,topChgTp_xz_C,actionTm_xz_C,topChgTm_xz_C
      double precision,dimension(0:ny-1,0:nz-1,0:nt-1,0:nloop,nL(that),offmax) :: actionT4p_xy_C,actionT4m_xy_C
      double precision,dimension(0:ny-1,0:nz-1,0:nt-1,0:nloop,nL(that),offmax) :: actionT4p_xz_C,actionT4m_xz_C
      double precision,dimension(0:ny-1,0:nz-1,0:nt-1,0:nloop,nL(that),offmax) :: topChgT4p_xy_C,topChgT4m_xy_C
      double precision,dimension(0:ny-1,0:nz-1,0:nt-1,0:nloop,nL(that),offmax) :: topChgT4p_xz_C,topChgT4m_xz_C
!HPF$ ALIGN (*,:,:,*,*,*) WITH correl(*,:,:,*,*) :: actionT4p_xy_C,topChgT4p_xy_C,actionT4m_xy_C,topChgT4m_xy_C
!HPF$ ALIGN (*,:,:,*,*,*) WITH correl(*,:,:,*,*) :: actionT4p_xz_C,topChgT4p_xz_C,actionT4m_xz_C,topChgT4m_xz_C

