c
c     Program to illustrate correlations between static quarks in a baryon and 
c     gluon action or topological charge
c     Deals with with T shaped Wilson loops corresponding to the Y shaped one.
c     Adapted from VacuumRespY_plan.f by F. Bissey, Oct. 2004
c
      program VacuumResp

      USE L_baryonParam
      USE GS_TADPOLEIMP
      USE GS_READLINKS
      USE L_ape_smear
      USE L_ape_smear3D
      USE L_fMuNu
      USE L_TopQandReconAction
      USE L_EleMagField
      USE L_EPSILONINDEX
      USE L_LOOPS
      USE L_WRITESHAPE

      IMPLICIT NONE

c     general variables
      integer,parameter                                          :: nsteps=4
      integer,parameter                                          :: offmax=4
      double precision, parameter                                :: alpha=0.7d0
c
      integer                                                    :: offset
      integer                                                    :: smear3d
      double precision                                           :: beta
      logical                                                    :: smear_links

!HPF$ TEMPLATE link(nx,ny,nz,nt,nc,nc)
!HPF$ DISTRIBUTE link(*,*,BLOCK,BLOCK,*,*)
!HPF$ TEMPLATE wloops(nx,ny,nz,nt,nL(that))
!HPF$ DISTRIBUTE wloops(BLOCK,BLOCK,*,*,*)
!HPF$ TEMPLATE correl(0:ny-1,0:nz-1,0:nt-1,nL(that),offmax)
!HPF$ DISTRIBUTE correl(*,BLOCK,BLOCK,*,*)

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
      double precision,dimension(nx,ny,nz,nt)                    :: reconAction
      double precision,dimension(nx,ny,nz,nt)                    :: topQDensity
!HPF$ ALIGN (:,:,*,*) WITH wloops(:,:,*,*,*) :: reconAction,topQDensity
      double precision,dimension(nx,ny,nz,nt)                    :: actionTx,actionTy,actionTz
!HPF$ ALIGN (:,:,*,*) WITH wloops(:,:,*,*,*) :: actionTx,actionTy,actionTz
      double precision,dimension(nL(that),offmax)                :: avgAction, avgTopChg
!HPF$ DISTRIBUTE avgAction(BLOCK,*)
!HPF$ DISTRIBUTE avgTopChg(BLOCK,*)
      double precision,dimension(nx,ny,nz,nt)                    :: actionT
      double precision,dimension(nx,ny,nz,nt)                    :: topChgT
!HPF$ ALIGN (:,:,*,*) WITH wloops(:,:,*,*,*) :: actionT,topChgT
      double precision,dimension(nx,ny,nz,nt)                    :: topChgTx,topChgTy,topChgTz
!HPF$ ALIGN (:,:,*,*) WITH wloops(:,:,*,*,*) :: topChgTx,topChgTy,topChgTz
c
c  Counting indexes
c
      integer                                                    :: ix,iy,iz,it,Tee,is
      integer                                                    :: icx,icy,icz
      integer                                                    :: incfg,ficfg,cfg
      integer                                                    :: ilp
c
c  Timer support
c
      integer                                                    :: start_count, end_count, count_rate
      integer                                                    :: time1, time2
      real                                                       :: elapsed_time

      character(len=132)                                         :: filename,configfile
c
c  loop dependent variables
c
c     if lattice s16t32
      integer, parameter                                      :: nloop=12 
c     if lattice s12t24
c      integer, parameter                                      :: nloop = 8
c
c  Time averaged Quantities
c
      double precision,dimension(nx,ny,nz,nt,0:nloop,nL(that)) :: WTp_xy,WTm_xy ! T-shape wilson loops
      double precision,dimension(nx,ny,nz,nt,0:nloop,nL(that)) :: WTp_xz,WTm_xz ! T-shape wilson loops
!HPF$ ALIGN (:,:,*,*,*,*,*) WITH wloops(:,:,*,*,*) :: WTp_xy,WTm_xy,WTp_xz,WTm_xz
      double precision,dimension(nx,ny,nz,nt,0:nloop,nL(that)) :: WT4p_xy,WT4m_xy ! T-shape wilson loops
      double precision,dimension(nx,ny,nz,nt,0:nloop,nL(that)) :: WT4p_xz,WT4m_xz ! T-shape wilson loops
!HPF$ ALIGN (:,:,*,*,*,*,*) WITH wloops(:,:,*,*,*) :: WT4p_xy,WT4m_xy,WT4p_xz,WT4m_xz

      double precision,dimension(0:nloop,nL(that))             :: WTp_xy_avg,WTm_xy_avg
      double precision,dimension(0:nloop,nL(that))             :: WTp_xz_avg,WTm_xz_avg
!HPF$ DISTRIBUTE WTp_xy_avg(BLOCK,*,*)
!HPF$ DISTRIBUTE WTm_xy_avg(BLOCK,*,*)
!HPF$ DISTRIBUTE WTp_xz_avg(BLOCK,*,*)
!HPF$ DISTRIBUTE WTm_xz_avg(BLOCK,*,*)
      double precision,dimension(0:nloop,nL(that))             :: WT4p_xy_avg,WT4m_xy_avg
      double precision,dimension(0:nloop,nL(that))             :: WT4p_xz_avg,WT4m_xz_avg
!HPF$ DISTRIBUTE WT4p_xy_avg(BLOCK,*,*)
!HPF$ DISTRIBUTE WT4m_xy_avg(BLOCK,*,*)
!HPF$ DISTRIBUTE WT4p_xz_avg(BLOCK,*,*)
!HPF$ DISTRIBUTE WT4m_xz_avg(BLOCK,*,*)
c
c  Correlations
c
      double precision,dimension(0:ny-1,0:nz-1,0:nt-1,0:nloop,nL(that),offmax) :: actionTp_xy_C,actionTm_xy_C
      double precision,dimension(0:ny-1,0:nz-1,0:nt-1,0:nloop,nL(that),offmax) :: actionTp_xz_C,actionTm_xz_C
      double precision,dimension(0:ny-1,0:nz-1,0:nt-1,0:nloop,nL(that),offmax) :: topChgTp_xy_C,topChgTm_xy_C
      double precision,dimension(0:ny-1,0:nz-1,0:nt-1,0:nloop,nL(that),offmax) :: topChgTp_xz_C,topChgTm_xz_C
!HPF$ ALIGN (*,:,:,*,*,*,*) WITH correl(*,:,:,*,*) :: actionTp_xy_C,topChgTp_xy_C,actionTm_xy_C,topChgTm_xy_C
!HPF$ ALIGN (*,:,:,*,*,*,*) WITH correl(*,:,:,*,*) :: actionTp_xz_C,topChgTp_xz_C,actionTm_xz_C,topChgTm_xz_C
      double precision,dimension(0:ny-1,0:nz-1,0:nt-1,0:nloop,nL(that),offmax) :: actionT4p_xy_C,actionT4m_xy_C
      double precision,dimension(0:ny-1,0:nz-1,0:nt-1,0:nloop,nL(that),offmax) :: actionT4p_xz_C,actionT4m_xz_C
      double precision,dimension(0:ny-1,0:nz-1,0:nt-1,0:nloop,nL(that),offmax) :: topChgT4p_xy_C,topChgT4m_xy_C
      double precision,dimension(0:ny-1,0:nz-1,0:nt-1,0:nloop,nL(that),offmax) :: topChgT4p_xz_C,topChgT4m_xz_C
!HPF$ ALIGN (*,:,:,*,*,*,*) WITH correl(*,:,:,*,*) :: actionT4p_xy_C,topChgT4p_xy_C,actionT4m_xy_C,topChgT4m_xy_C
!HPF$ ALIGN (*,:,:,*,*,*,*) WITH correl(*,:,:,*,*) :: actionT4p_xz_C,topChgT4p_xz_C,actionT4m_xz_C,topChgT4m_xz_C

c
c  Begin execution
c
      write(configfile,'(a)') "./dqconfig"


      pi = 4.0d0 * atan(1.0d0)
      volume = nx * ny * nz * nt

      open(101,file=configfile)
      read (101,'(a132)') basedir
      write(*,'(a)') basedir(1:len_trim(basedir))

      read (101,'(a132)') lastConfig
      write(*,'(a)') lastConfig(1:len_trim(lastConfig))

      read (101,*) incfg
      write(*,*) incfg

      read (101,*) ficfg
      write(*,*) ficfg

      read (101,'(a132)') prefix
      write(*,'(a)') prefix(1:len_trim(prefix))

      read (101,*) Qpaths
      write(*,*) Qpaths

      read (101,*) smear3d
      write(*,*) smear3d

      read (101,*) Correlate
      write(*,*) Correlate

      read(101,'(l7)') smear_links
      write(*,'(l7)') smear_links
      write(*,*)

      close(101)
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     initialize the epsilon index for the calls from lt_oop
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      call setEpsilonIndex()

      do cfg = incfg, ficfg

         CALL SYSTEM_CLOCK(start_count, count_rate)

         write(input,fmt='(a,a,i3.3)')trim(basedir),trim(lastconfig),cfg

         call ReadLinks(input,ur,ui,nfig,beta,nxold,nyold,nzold,ntold,lastPlaq,plaqbarAvg,uzero)

         if (nx.ne.nxold) then
            write(*,'(a)') 'Mismatch in nx.'
            goto 999
         else if (ny.ne.nyold) then
            write(*,'(a)') 'Mismatch in ny.'
            goto 999
         else if (nz.ne.nzold) then
            write(*,'(a)') 'Mismatch in nz.'
            goto 999
         else if (nt.ne.ntold) then
            write(*,'(a)') 'Mismatch in nt.'
            goto 999
         else if ( nfig.eq.0 ) then
            write(*,'(a)') 'Mismatch in average plaquette.'
            go to 999
         end if

         if ( smear_links ) then

            write (output,'(a,a,i3.3,a,i2.2,a,i2.2)')trim(prefix),trim(lastconfig),cfg,'.s',nsteps,'-s',smear3d

         else

            write(output,'(a,a,i3.3,a,i2.2)')trim(prefix),trim(lastconfig),cfg,'.s00-s',smear3d

         endif

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  The gauge field has been successfully read
c  Get the field strength tensor
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

         utr = ur
         uti = ui
         if ( smear_links ) then
            do istep = 1,nsteps

               call ape_smear(utr,uti,alpha)

            end do
         end if

         call tadpoleimp(utr,uti,uzero)
         call fMuNu(utr,uti,Fr,Fi,uzero,Qpaths)

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

         if (Correlate == 1) then

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  Calculation of topological charge density
c  Scale results to single instanton results
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

            call TopQandReconAction(Fr,Fi,uzero,topQDensity,Q,reconAction,pluckbar)

            sons1 = beta * pluckbar * nx * ny * nz * nt * (mu*(mu-1)/2.0d0) /
     &           ( 8.0d0*pi**2/(6.0d0/beta) )

            reconAction = reconAction * beta * (mu*(mu-1)/2.0d0) / ( 8.0d0*pi**2/(6.0d0/beta) )

            write(*,'(/,2(a,f15.8,/),/)')
     &           'Reconstructed S/S_0 = ', sons1,
     &           'Topological Charge  = ', Q

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

         else

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  We'll put the electric field in the reconAction and the
c                magnetic field in the topQdensity
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

            call EleMagField(Fr,Fi,reconAction,topQdensity)

            write(*,'(/,2(a,e15.8,/),/)')
     &           'Sum of Electric Field = ', SUM(reconAction),
     &           'Sum of Magnetic Field = ', SUM(topQdensity)

         end if

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  Get the Wilson Loops
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

         utr = ur
         uti = ui

         do istep = 1,smear3d

            call ape_smear3D(utr,uti,alpha,that)

         end do

         if ( Correlate == 1) then

            fstr1='.action.correl.unf'
            fstr2='.topChg.correl.unf'

         else

            fstr1='.ele.correl.unf'
            fstr2='.mag.correl.unf'

         endif

         call system_clock(time1)

      call t_loops(xhat,yhat,utr,uti,
     &             WTp_xy(:,:,:,:,:,:),WTm_xy(:,:,:,:,:,:),WT4p_xy(:,:,:,:,:,:),WT4m_xy(:,:,:,:,:,:))
      call t_loops(xhat,zhat,utr,uti,
     &             WTp_xz(:,:,:,:,:,:),WTm_xz(:,:,:,:,:,:),WT4p_xz(:,:,:,:,:,:),WT4m_xz(:,:,:,:,:,:))

      do ilp = 1, nloop

         do it = 1, nL(that)

            WTp_xy_avg(ilp,it)  = Sum( WTp_xy(:,:,:,:,ilp,it) )  / volume
            WTm_xy_avg(ilp,it)  = Sum( WTm_xy(:,:,:,:,ilp,it) )  / volume
            WTp_xz_avg(ilp,it)  = Sum( WTp_xz(:,:,:,:,ilp,it) )  / volume
            WTm_xz_avg(ilp,it)  = Sum( WTm_xz(:,:,:,:,ilp,it) )  / volume
            WT4p_xy_avg(ilp,it) = Sum( WT4p_xy(:,:,:,:,ilp,it) ) / volume
            WT4m_xy_avg(ilp,it) = Sum( WT4m_xy(:,:,:,:,ilp,it) ) / volume
            WT4p_xz_avg(ilp,it) = Sum( WT4p_xz(:,:,:,:,ilp,it) ) / volume
            WT4m_xz_avg(ilp,it) = Sum( WT4m_xz(:,:,:,:,ilp,it) ) / volume

         end do              ! it = 2, nL(that)  loop end
      end do

         call system_clock(time2)
         write(*,'(a,f15.8,a)') 'time spent constructing loops is: ',(real(time2-time1)/real(count_rate)),' s'
         write(*,*) offmax
c
         call system_clock(time1)
         do offset = 1,offmax
            do Tee = 2*offset,nL(that)
               do it = 1,nx

                  actionT(it,:,:,:) = 0.0d0
                  topChgT(it,:,:,:) = 0.0d0

                  do is = it+offset , it+Tee-offset

                     actionT(it,:,:,:) = actionT(it,:,:,:) + reconaction(( is - ((is - 1)/nx)*nx ),:,:,:)

                     topChgT(it,:,:,:) = topChgT(it,:,:,:) + abs( topQDensity(( is - ((is - 1)/nx)*nx ),:,:,:) )

                  end do           ! is loop end 

                  actionT(it,:,:,:) = actionT(it,:,:,:)/(Tee - 2*offset + 1)
                  topChgT(it,:,:,:) = topChgT(it,:,:,:)/(Tee - 2*offset + 1)

               end do              ! it = 1, nx loop end

               avgAction(Tee,offset) = Sum( actionT(:,:,:,:) ) / volume
               avgTopChg(Tee,offset) = Sum( topChgT(:,:,:,:) ) / volume

               do iy = 0, ny-1

                  actionTx(:,:,:,:) = cshift(actionT(:,:,:,:), dim = 2, shift = iy)
                  topChgTx(:,:,:,:) = cshift(topChgT(:,:,:,:), dim = 2, shift = iy)

                  do iz = 0, nz-1

                     actionTy(:,:,:,:) = cshift(actionTx(:,:,:,:), dim = 3, shift = iz)
                     topChgTy(:,:,:,:) = cshift(topChgTx(:,:,:,:), dim = 3, shift = iz)

                     do it = 0, nt-1

                        actionTz(:,:,:,:) = cshift(actionTy(:,:,:,:), dim = 4, shift = it)
                        topChgTz(:,:,:,:) = cshift(topChgTy(:,:,:,:), dim = 4, shift = it)

      do ilp = 1, nloop

         actionTp_xy_C(iy,iz,it,ilp,Tee,offset) = sum( WTp_xy(:,:,:,:,ilp,Tee) * actionTz(:,:,:,:) ) / volume
         actionTm_xy_C(iy,iz,it,ilp,Tee,offset) = sum( WTm_xy(:,:,:,:,ilp,Tee) * actionTz(:,:,:,:) ) / volume
         actionTp_xz_C(iy,iz,it,ilp,Tee,offset) = sum( WTp_xz(:,:,:,:,ilp,Tee) * actionTz(:,:,:,:) ) / volume
         actionTm_xz_C(iy,iz,it,ilp,Tee,offset) = sum( WTm_xz(:,:,:,:,ilp,Tee) * actionTz(:,:,:,:) ) / volume

         topChgTp_xy_C(iy,iz,it,ilp,Tee,offset) = sum( WTp_xy(:,:,:,:,ilp,Tee) * topChgTz(:,:,:,:) ) / volume
         topChgTm_xy_C(iy,iz,it,ilp,Tee,offset) = sum( WTm_xy(:,:,:,:,ilp,Tee) * topChgTz(:,:,:,:) ) / volume
         topChgTp_xz_C(iy,iz,it,ilp,Tee,offset) = sum( WTp_xz(:,:,:,:,ilp,Tee) * topChgTz(:,:,:,:) ) / volume
         topChgTm_xz_C(iy,iz,it,ilp,Tee,offset) = sum( WTm_xz(:,:,:,:,ilp,Tee) * topChgTz(:,:,:,:) ) / volume

         actionT4p_xy_C(iy,iz,it,ilp,Tee,offset) = sum( WT4p_xy(:,:,:,:,ilp,Tee) * actionTz(:,:,:,:) ) / volume
         actionT4m_xy_C(iy,iz,it,ilp,Tee,offset) = sum( WT4m_xy(:,:,:,:,ilp,Tee) * actionTz(:,:,:,:) ) / volume
         actionT4p_xz_C(iy,iz,it,ilp,Tee,offset) = sum( WT4p_xz(:,:,:,:,ilp,Tee) * actionTz(:,:,:,:) ) / volume
         actionT4m_xz_C(iy,iz,it,ilp,Tee,offset) = sum( WT4m_xz(:,:,:,:,ilp,Tee) * actionTz(:,:,:,:) ) / volume

         topChgT4p_xy_C(iy,iz,it,ilp,Tee,offset) = sum( WT4p_xy(:,:,:,:,ilp,Tee) * topChgTz(:,:,:,:) ) / volume
         topChgT4m_xy_C(iy,iz,it,ilp,Tee,offset) = sum( WT4m_xy(:,:,:,:,ilp,Tee) * topChgTz(:,:,:,:) ) / volume
         topChgT4p_xz_C(iy,iz,it,ilp,Tee,offset) = sum( WT4p_xz(:,:,:,:,ilp,Tee) * topChgTz(:,:,:,:) ) / volume
         topChgT4m_xz_C(iy,iz,it,ilp,Tee,offset) = sum( WT4m_xz(:,:,:,:,ilp,Tee) * topChgTz(:,:,:,:) ) / volume

      end do     ! ilp loop end

                     end do        ! it = 0, nt-1 loop end
                  end do           ! iz loop end
               end do              ! iy loop end
            end do                 ! Tee loop end
         end do                    ! offset loop end
c
         call system_clock(time2)


      do ix = 0, nt-1

         icx = modulo(nt - ix , nt) 

         actionTp_xy_C(:,:,ix,1:nloop,:,:) = actionTp_xy_C(:,:, ix,1:nloop,:,:) + actionTm_xy_C(:,:,icx,1:nloop,:,:)

         topChgTp_xy_C(:,:,ix,1:nloop,:,:) = topChgTp_xy_C(:,:, ix,1:nloop,:,:) + topChgTm_xy_C(:,:,icx,1:nloop,:,:)

         actionT4p_xy_C(:,:,ix,1:nloop,:,:) = actionT4p_xy_C(:,:, ix,1:nloop,:,:) + actionT4m_xy_C(:,:,icx,1:nloop,:,:)

         topChgT4p_xy_C(:,:,ix,1:nloop,:,:) = topChgT4p_xy_C(:,:, ix,1:nloop,:,:) + topChgT4m_xy_C(:,:,icx,1:nloop,:,:)

      end do
c
      WTp_xy_avg(1:nloop,:) = WTp_xy_avg(1:nloop,:) + WTm_xy_avg(1:nloop,:)
      WT4p_xy_avg(1:nloop,:) = WT4p_xy_avg(1:nloop,:) + WT4m_xy_avg(1:nloop,:)

      do iy = 0, nz-1
         do iz = 0, ny-1

            icy = modulo(nz - iz , nz)
            icz = modulo(     iy , ny)

            actionTp_xy_C(iz,iy,:,1:nloop,:,:) = actionTp_xy_C(iz,iy,:,1:nloop,:,:) + actionTp_xz_C(icz,icy,:,1:nloop,:,:)

            topChgTp_xy_C(iz,iy,:,1:nloop,:,:) = topChgTp_xy_C(iz,iy,:,1:nloop,:,:) + topChgTp_xz_C(icz,icy,:,1:nloop,:,:)

            actionT4p_xy_C(iz,iy,:,1:nloop,:,:) = actionT4p_xy_C(iz,iy,:,1:nloop,:,:) + actionT4p_xz_C(icz,icy,:,1:nloop,:,:)

            topChgT4p_xy_C(iz,iy,:,1:nloop,:,:) = topChgT4p_xy_C(iz,iy,:,1:nloop,:,:) + topChgT4p_xz_C(icz,icy,:,1:nloop,:,:)

         end do
      end do
c
      WTp_xy_avg(1:nloop,:) = WTp_xy_avg(1:nloop,:) + WTp_xz_avg(1:nloop,:)
      WT4p_xy_avg(1:nloop,:) = WT4p_xy_avg(1:nloop,:) + WT4p_xz_avg(1:nloop,:)

      do iy = 0, nz-1
         do iz = 0, ny-1
            do ix = 0, nt-1

               icx = modulo(nt - ix , nt) 
               icy = modulo(nz - iz , nz)
               icz = modulo(     iy , ny)

               actionTp_xy_C(iz,iy,ix,1:nloop,:,:) = actionTp_xy_C(iz,iy,ix,1:nloop,:,:) + actionTm_xz_C(icz,icy,icx,1:nloop,:,:)

               topChgTp_xy_C(iz,iy,ix,1:nloop,:,:) = topChgTp_xy_C(iz,iy,ix,1:nloop,:,:) + topChgTm_xz_C(icz,icy,icx,1:nloop,:,:)

               actionT4p_xy_C(iz,iy,ix,1:nloop,:,:) = actionT4p_xy_C(iz,iy,ix,1:nloop,:,:) + actionT4m_xz_C(icz,icy,icx,1:nloop,:,:)

               topChgT4p_xy_C(iz,iy,ix,1:nloop,:,:) = topChgT4p_xy_C(iz,iy,ix,1:nloop,:,:) + topChgT4m_xz_C(icz,icy,icx,1:nloop,:,:)

            end do
         end do
      end do
c
      WTp_xy_avg(1:nloop,:) = WTp_xy_avg(1:nloop,:) + WTm_xz_avg(1:nloop,:)
      WT4p_xy_avg(1:nloop,:) = WT4p_xy_avg(1:nloop,:) + WT4m_xz_avg(1:nloop,:)

      actionTp_xy_C(:,:,:,1:nloop,:,:) = actionTp_xy_C(:,:,:,1:nloop,:,:) / 4
      topChgTp_xy_C(:,:,:,1:nloop,:,:) = topChgTp_xy_C(:,:,:,1:nloop,:,:) / 4
      WTp_xy_avg(1:nloop,:) = WTp_xy_avg(1:nloop,:) / 4

      actionT4p_xy_C(:,:,:,1:nloop,:,:) = actionT4p_xy_C(:,:,:,1:nloop,:,:) / 4
      topChgT4p_xy_C(:,:,:,1:nloop,:,:) = topChgT4p_xy_C(:,:,:,1:nloop,:,:) / 4
      WT4p_xy_avg(1:nloop,:) = WT4p_xy_avg(1:nloop,:) / 4

      actionTp_xy_C = cshift(cshift(actionTp_xy_C,dim=3,shift= 4),dim=2,shift= ny/2)
      actionT4p_xy_C = cshift(cshift(actionT4p_xy_C,dim=3,shift= 4),dim=2,shift= ny/2)

      topChgTp_xy_C = cshift(cshift(topChgTp_xy_C,dim=3,shift= 4),dim=2,shift= ny/2)
      topChgT4p_xy_C = cshift(cshift(topChgT4p_xy_C,dim=3,shift= 4),dim=2,shift= ny/2)
      filename=trim(output)//'.DQ2avg4'//trim(fstr1)
      call writeDQshape(filename,actionTp_xy_C(:,:,:,:,:,:),WTp_xy_avg(:,:),avgAction,offmax)

      filename=trim(output)//'.DQ2avg4'//trim(fstr2)
      call writeDQshape(filename,TopChgTp_xy_C(:,:,:,:,:,:),WTp_xy_avg(:,:),avgTopChg,offmax)

      filename=trim(output)//'.DQ4avg4'//trim(fstr1)
      call writeDQshape(filename,actionT4p_xy_C(:,:,:,:,:,:),WT4p_xy_avg(:,:),avgAction,offmax)

      filename=trim(output)//'.DQ4avg4'//trim(fstr2)
      call writeDQshape(filename,TopChgT4p_xy_C(:,:,:,:,:,:),WT4p_xy_avg(:,:),avgTopChg,offmax)

 999     CALL SYSTEM_CLOCK(end_count)

         elapsed_time = REAL(end_count - start_count) / REAL(count_rate)
c
c  Excluding assignmentgs
c
         if (elapsed_time .ne. 0.0) then

            write(*,'(/,a,f20.5,a,i3.3)') 'The elapsed time is',elapsed_time,
     &           ' seconds in configuration:',cfg

         endif

      end do                    !end cfg loop

      end program VacuumResp

