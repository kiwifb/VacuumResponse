c     general variables
      integer,parameter                                          :: offmax=4
      integer,parameter                                          :: nsteps=4
      double precision, parameter                                :: alpha=0.7d0
c
      integer                                                    :: offset
      integer                                                    :: smear3d
      double precision                                           :: beta     
      logical                                                    :: smear_links
!HPF$ TEMPLATE link(nx,ny,nz,nt,nc,nc)
!HPF$ DISTRIBUTE link(*,*,BLOCK,BLOCK,*,*)
      double precision,dimension(nx,ny,nz,nt,mu,nc,nc)           :: ur,ui,utr,uti
!HPF$ ALIGN (*,*,:,:,*,*,*) WITH link(*,*,:,:,*,*) :: ur,ui,utr,uti
c
c     local variables
c
      integer,parameter                                          :: nf=6         !# of Fmunu
      double precision                                           :: pi
      integer                                                    :: Qpaths,Correlate,volume
      integer                                                    :: nfig,istep
      double precision                                           :: uzero=1.0d0
      double precision                                           :: pluckbar
      double precision                                           :: sons1
      double precision                                           :: plaqbarAvg,lastPlaq
      character(len=132)                                         :: lastconfig,prefix,basedir
      character(len=132)                                         :: input,output
      character(len=50)                                          :: fstr1,fstr2
      integer                                                    :: nxold,nyold,nzold,ntold
      double precision,dimension(nx,ny,nz,nt,nc,nc,nf)           :: Fr,Fi        !Fmunu
!HPF$ ALIGN (*,*,:,:,*,*,*) WITH link(*,*,:,:,*,*) :: Fr,Fi
      double precision                                           :: Q=0.0d0 
!HPF$ TEMPLATE wloops(nx,ny,nz,nt,nL(that))
!HPF$ DISTRIBUTE wloops(BLOCK,BLOCK,*,*,*)
      double precision,dimension(nx,ny,nz,nt)                    :: reconAction
      double precision,dimension(nx,ny,nz,nt)                    :: topQDensity
!HPF$ ALIGN (:,:,*,*) WITH wloops(:,:,*,*,*) :: reconAction,topQDensity
      double precision,dimension(nx,ny,nz,nt)                    :: actionTx,actionTy,actionTz
!HPF$ ALIGN (:,:,*,*) WITH wloops(:,:,*,*,*) :: actionTx,actionTy,actionTz
      double precision,dimension(nL(that),offmax)                :: avgAction, avgTopChg
!HPF$ DISTRIBUTE avgAction(BLOCK,*)
!HPF$ DISTRIBUTE avgTopChg(BLOCK,*)
      double precision,dimension(nx,ny,nz,nt)                    :: topChgTx,topChgTy,topChgTz
!HPF$ ALIGN (:,:,*,*) WITH wloops(:,:,*,*,*) :: topChgTx,topChgTy,topChgTz
!HPF$ TEMPLATE correl(0:ny-1,0:nz-1,0:nt-1,nL(that),offmax)
!HPF$ DISTRIBUTE correl(*,BLOCK,BLOCK,*,*)
c
c  Counting indexes
c
      integer                                                    :: iy,iz,it,Tee,is
      integer                                                    :: incfg,ficfg,cfg
c
c  Timer support
c
      integer                                                    :: start_count, end_count, count_rate
      integer                                                    :: time1, time2
      real                                                       :: elapsed_time
