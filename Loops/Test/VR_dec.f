c
c  Time averaged Quantities
c
      double precision,dimension(nx,ny,nz,nt,0:nloop,nL(that))  :: WTp ! T-shape wilson loops
!HPF$ ALIGN (:,:,*,*,*,*) WITH wloops(:,:,*,*,*) :: WTp

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Nomenclature for the T shape Wilson loop
c     The index running from 0 to nY correspond to the various
c     position of the three quarks in the T-loop.
c     The sets of coordinates for each value of the index are in ty-loops.f
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      double precision,dimension(0:nloop,nL(that))            :: WTpavg
!HPF$ DISTRIBUTE WTpavg(BLOCK,*)
c
c  Correlations
c
      double precision,dimension(0:ny-1,0:nz-1,0:nt-1,0:nloop,nL(that),offmax) :: actionTpC
!HPF$ ALIGN (*,:,:,*,*,*) WITH correl(*,:,:,*,*) :: actionTpC
