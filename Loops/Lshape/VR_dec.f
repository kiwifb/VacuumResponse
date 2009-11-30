c
c  Time averaged Quantities
c
      double precision,dimension(nx,ny,nz,nt,0:nloop,4,nL(that))   :: WLxy,WLxz ! T-shape wilson loops
!HPF$ ALIGN (:,:,*,*,*,*,*) WITH Wloops(:,:,*,*,*) :: WLxy,WLxz

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Nomenclature for the L shape Wilson loop
c     The index running from 0 to nY correspond to the various
c     position of the three quarks in the L-loop.
c     The loop #1 in the first quadrant is made of quarks in the
c     following positions: (0,0), (0,i) and (i,0)
c     #2 (0,0), (0,i) and (-i,0)
c     #3 (0.0), (0,-i) and (-i,0)
c     #4 (0,0), (0,-i) and (i,0)
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      double precision,dimension(0:nloop,nL(that))               :: WLavg
!HPF$ DISTRIBUTE WLavg(BLOCK,*)
      double precision,dimension(nx,ny,nz,nt,nc,nc)              :: uptr,upti
!HPF$ DISTRIBUTE uptr(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE upti(*,*,BLOCK,BLOCK,*,*)
c
c  Correlations
c
      double precision,dimension(0:ny-1,0:nz-1,0:nt-1,0:nloop,nL(that),offmax) :: actionLC
      double precision,dimension(0:ny-1,0:nz-1,0:nt-1,0:nloop,nL(that),offmax) :: topChgLC
!HPF$ DISTRIBUTE actionLC(*,*,*,BLOCK,BLOCK,*)
!HPF$ DISTRIBUTE topChgLC(*,*,*,BLOCK,BLOCK,*)
