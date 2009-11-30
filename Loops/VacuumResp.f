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
      USE L_product

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

      include'Templates/Wloop_template.f'

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
      integer                                                    :: icx,icy,icz,ict
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
#include"loopsize.f"
#include"VR_dec.f"
c
c  Begin execution
c
#include"VR_config.f"

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

#include"VR_loop.f"

         call system_clock(time2)
         write(*,'(a,f15.8,a)') 'time spent constructing loops is: ',(real(time2-time1)/real(count_rate)),' s'
         write(*,*) offmax
c
         call system_clock(time1)
!$omp parallel do
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

#include"VR_correl.f"

                     end do        ! it = 0, nt-1 loop end
                  end do           ! iz loop end
               end do              ! iy loop end
            end do                 ! Tee loop end
         end do                    ! offset loop end
!$omp end parallel do
c
         call system_clock(time2)

#include"VR_write.f"

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

