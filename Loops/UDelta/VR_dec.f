c
c  Time averaged Quantities
c
      double precision,dimension(nx,ny,nz,nt,0:nloop,nL(that))  :: WTp_xy, WTp_xz! T-shape wilson loops
!HPF$ ALIGN (:,:,*,*,*,*) WITH wloops(:,:,*,*,*) :: WTpxy
!HPF$ ALIGN (:,:,*,*,*,*) WITH wloops(:,:,*,*,*) :: WTpxz

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Nomenclature for the T shape Wilson loop
c     The index running from 0 to nY correspond to the various
c     position of the three quarks in the T-loop.
c     The sets of coordinates for each value of the index are in ty-loops.f
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      double precision,dimension(0:nloop,nL(that))            :: WTpavg_xy, WTpavg_xz
!HPF$ DISTRIBUTE WTpavg_xy(BLOCK,*)
!HPF$ DISTRIBUTE WTpavg_xz(BLOCK,*)
c
c  Correlations
c
      double precision,dimension(0:ny-1,0:nz-1,0:nt-1,0:nloop,nL(that),offmax) :: actionTpC_xy, actionTpC_xz
      double precision,dimension(0:ny-1,0:nz-1,0:nt-1,0:nloop,nL(that),offmax) :: topChgTpC_xy, topChgTpC_xz
!HPF$ ALIGN (*,:,:,*,*,*) WITH correl(*,:,:,*,*) :: actionTpC_xy,topChgTpC_xy
!HPF$ ALIGN (*,:,:,*,*,*) WITH correl(*,:,:,*,*) :: actionTpC_xz,topChgTpC_xz
