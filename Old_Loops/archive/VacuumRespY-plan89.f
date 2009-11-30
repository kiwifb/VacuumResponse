c
c     Program to illustrate correlations between static quarks in a baryon and 
c     gluon action or topological charge
c     Deals with with Y shaped Wilson loops.
c     Adapted from VacuumRespLT.f by F. Bissey, Jan. 2004
c
      program VacuumRespY

      USE L_baryonParam
      USE GS_TADPOLEIMP
      USE GS_READLINKS
      USE L_ape_smear
      USE L_ape_smear3D
      USE L_fMuNu
      USE L_TopQandReconAction
      USE L_EleMagField
      USE L_YLOOPS89
      USE L_EPSILONINDEX
      USE L_WRITEYSHAPE89

      IMPLICIT NONE

c     general variables
      integer,parameter                                          :: offmax=4
      integer,parameter                                          :: nsteps=4
      double precision, parameter                                :: alpha=0.7d0
c
      integer                                                    :: offset
      double precision                                           :: beta     
      logical                                                    :: smear_links
      double precision,dimension(nx,ny,nz,nt,mu,nc,nc)           :: ur,ui,utr,uti
!HPF$ DISTRIBUTE ur(*,*,BLOCK,BLOCK,*,*,*)
!HPF$ ALIGN (*,*,:,:,*,*,*) WITH ur(*,*,:,:,*,*,*) :: ui, utr, uti
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
      character(len=132)                                         :: lastconfig,string,prefix,basedir
      character(len=132)                                         :: input,output
      character(len=50)                                          :: fstr1,fstr2
      character(len=132)                                         :: filename
      integer                                                    :: nxold,nyold,nzold,ntold
      double precision,dimension(nx,ny,nz,nt,nc,nc,nf)           :: Fr,Fi        !Fmunu
!HPF$ DISTRIBUTE Fr(*,*,BLOCK,BLOCK,*,*,*)
!HPF$ DISTRIBUTE Fi(*,*,BLOCK,BLOCK,*,*,*)
      double precision,dimension(nx,ny,nz,nt)                    :: reconAction
!HPF$ DISTRIBUTE reconAction(BLOCK,BLOCK,*,*)
      double precision,dimension(nx,ny,nz,nt)                    :: topQDensity
!HPF$ ALIGN (:,:,*,*) WITH reconAction(:,:,*,*) :: topQDensity
      double precision                                           :: Q=0.0d0 
c
c  Time averaged Quantities
c
      double precision,dimension(nx,ny,nz,nt,8:nYloop,nL(that),2)  :: WYp,WYm ! Y-shape wilson loops
!HPF$ ALIGN (:,:,*,*,*,*,*) WITH reconAction(:,:,*,*) :: WYp,WYm

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Nomenclature for the Y shape Wilson loop
c     The index running from 0 to nY correspond to the various
c     position of the three quarks in the Y-loop.
c     The sets of coordinates are in the file y-loops
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      double precision,dimension(8:nYloop,nL(that),2)            :: WYpavg,WYmavg
!HPF$ DISTRIBUTE WYpavg(*,BLOCK,*)
!HPF$ DISTRIBUTE WYmavg(*,BLOCK,*)
      double precision,dimension(nx,ny,nz,nt)                    :: actionT
      double precision,dimension(nx,ny,nz,nt)                    :: actionTx,actionTy,actionTz
!HPF$ ALIGN (:,:,*,*) WITH WYp(:,:,*,*,*,*,*) :: actionT,actionTx,actionTy,actionTz
      double precision,dimension(nL(that),offmax)                :: avgAction, avgTopChg
!HPF$ DISTRIBUTE avgAction(BLOCK,*)
!HPF$ DISTRIBUTE avgTopChg(BLOCK,*)
      double precision,dimension(nx,ny,nz,nt)                    :: topChgT
      double precision,dimension(nx,ny,nz,nt)                    :: topChgTx,topChgTy,topChgTz
!HPF$ ALIGN (:,:,*,*) WITH WYp(:,:,*,*,*,*,*) :: topChgT,topChgTx,topChgTy,topChgTz
c
c  Correlations
c
      double precision,dimension(0:ny-1,0:nz-1,0:nt-1,8:nYloop,nL(that),offmax,2) :: actionYpC,actionYmC
!HPF$ DISTRIBUTE actionYpC(*,BLOCK,BLOCK,*,*,*,*)
!HPF$ ALIGN (*,:,:,*,*,*,*) WITH actionYpC(*,:,:,*,*,*,*) :: actionYmC
      double precision,dimension(0:ny-1,0:nz-1,0:nt-1,8:nYloop,nL(that),offmax,2) :: topChgYpC,topChgYmc
!HPF$ ALIGN (*,:,:,*,*,*,*) WITH actionYpC(*,:,:,*,*,*,*) :: topChgYpC,topChgYmc
c
c  Counting indexes
c
      integer                                                    :: iy,iz,it,ic,Tee,is,iylp,iydir
      integer                                                    :: incfg,ficfg,cfg
      integer,dimension(2)                                       :: yplane
c
c  Timer support
c
      integer                                                    :: start_count, end_count, count_rate
      integer                                                    :: time1, time2
      real                                                       :: elapsed_time

c
c  Begin execution
c
      pi = 4.0d0 * atan(1.0d0)
      volume = nx * ny * nz * nt

      write(*,*)
      write(*,'(a)')'Please enter the configurations directory.'
      read (*,'(a132)') basedir
      write(*,'(a)') basedir(1:len_trim(basedir))

      write(*,*)
      write(*,'(a)')'Please enter the gauge field file name.'
      read (*,'(a132)') lastConfig
      write(*,'(a)') lastConfig(1:len_trim(lastConfig))

      write(*,*)
      write(*,'(a)')'Please enter the index of the first configuration.'
      read (*,*) incfg
      write(*,*) incfg

      write(*,*)
      write(*,'(a)')'Please enter the index of the last configuration.'
      read (*,*) ficfg
      write(*,*) ficfg

      write(*,*)
      write(*,'(a)')'Please enter directory to output to. (baryonCorrel/) .'
      write(*,'(a)')'Note: this subdirectory must exist.'
      read (*,'(a132)') prefix
      write(*,'(a)') prefix(1:len_trim(prefix))

      write(*,*)
      write(*,'(a)')'Which F_mu_nu do you desire? (3)'
      write(*,'(a)')'      1: Standard plaquette F_mu_nu'
      write(*,'(a)')'      2: (1x1), (1x2) F_mu_nu'
      write(*,'(a)')'      3: (1x1), (2x2), (3x3) F_mu_nu'
      write(*,'(a)')'      4: (1x1), (2x2), (1x2),(1X3) F_mu_nu (4LTQ #1)'
      write(*,'(a)')'      5: 5LI F_mu_nu'
      read (*,*) Qpaths
      write(*,*) Qpaths

      write(*,*)
      write(*,'(a)')'What would you like to correlate? (1)'
      write(*,'(a)')'      1: Action and topological charge.'
      write(*,'(a)')'      2: Electric and Magnetic fields.'
      read (*,*) Correlate
      write(*,*) Correlate
      
      write(*,*)'Smear Links For The FMuNu Subroutine? (.true. or .false.) (.true.)'
      read(*,'(l7)') smear_links
      write(*,'(l7)') smear_links
      write(*,*)
c'
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     initialize the epsilon index for the calls from y_loop
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      call setEpsilonIndex()        

      yplane(1)=yhat
      yplane(2)=zhat

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

            write (output,'(a,a,i3.3,a,i2.2)')trim(prefix),trim(lastconfig),cfg,'.s',nsteps

         else

            write(output,'(a,a,i3.3,a)')trim(prefix),trim(lastconfig),cfg,'.s00'

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

            call EleMagField(Fr,Fi,uzero,reconAction,topQdensity)

            write(*,'(/,2(a,e15.8,/),/)') 
     &           'Sum of Electric Field = ', SUM(reconAction), 
     &           'Sum of Magnetic Field = ', SUM(topQdensity)

         end if

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  Get the Wilson Loops
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

         utr = ur
         uti = ui

         do istep = 1,10
            
            call ape_smear3D(utr,uti,alpha,that)

         end do


         if ( Correlate == 1) then
               
            fstr1='89.action.correl.unf'
            fstr2='89.topChg.correl.unf'
            
         else

            fstr1='89.ele.correl.unf'
            fstr2='89.mag.correl.unf'
            
         endif

         xdir = xhat

         do iydir = 1,2
            
            ydir = yplane(iydir)

            call system_clock(time1)
            call y_loops_mirror(utr,uti,WYp(:,:,:,:,:,:,iydir),WYm(:,:,:,:,:,:,iydir))

            do iylp = 8, nYloop

               do it = 2, nL(that)
               
                  WYpavg(iylp,it,iydir) = Sum( WYp(:,:,:,:,iylp,it,iydir) ) / volume
                  WYmavg(iylp,it,iydir) = Sum( WYm(:,:,:,:,iylp,it,iydir) ) / volume

               end do           ! it = 2, nL(that)  loop end
            end do
         end do                 ! iydir loop end

         call system_clock(time2)
         write(*,'(a,f15.8,a)') 'time spent constructing Y loops is: ',(real(time2-time1)/real(count_rate)),' s'
         write(*,*) offmax
c     '
         call system_clock(time1)
         do offset = 1,offmax
            do Tee = 2*offset,nL(that)
               do it = 1,nx

                  actionT(it,:,:,:) = 0.0d0
                  topChgT(it,:,:,:) = 0.0d0

                  do is = it+offset , it+Tee-offset
                     
                     actionT(it,:,:,:) = actionT(it,:,:,:) + reconaction(( is - ((is - 1)/nx)*nx ),:,:,:)
                  
                     topChgT(it,:,:,:) = topChgT(it,:,:,:) + abs( topQDensity(( is - ((is - 1)/nx)*nx ),:,:,:) )

                  end do        ! is loop end 

                  actionT(it,:,:,:) = actionT(it,:,:,:)/(Tee - 2*offset + 1)
                  topChgT(it,:,:,:) = topChgT(it,:,:,:)/(Tee - 2*offset + 1)

               end do           ! it = 1, nx loop end

               avgAction(Tee,offset) = Sum( actionT(:,:,:,:) ) / volume
               avgTopChg(Tee,offset) = Sum( topChgT(:,:,:,:) ) / volume

               do iy = 0, ny-1

c     write(*,'(a,i2,a,i2)') 'Step ',iy,' of ',ny-1
                     
                  actionTx(:,:,:,:) = cshift(actionT(:,:,:,:), dim = 2, shift = iy)
                  topChgTx(:,:,:,:) = cshift(topChgT(:,:,:,:), dim = 2, shift = iy)
                  
                  do iz = 0, nz-1

                     actionTy(:,:,:,:) = cshift(actionTx(:,:,:,:), dim = 3, shift = iz)
                     topChgTy(:,:,:,:) = cshift(topChgTx(:,:,:,:), dim = 3, shift = iz)

                     do it = 0, nt-1
                     
                        actionTz(:,:,:,:) = cshift(actionTy(:,:,:,:), dim = 4, shift = it)
                        topChgTz(:,:,:,:) = cshift(topChgTy(:,:,:,:), dim = 4, shift = it)

                        do iylp = 8, nYloop

                           actionYpC(iy,iz,it,iylp,Tee,offset,1) = sum( WYp(:,:,:,:,iylp,Tee,1) * actionTz(:,:,:,:) ) / volume
                           actionYmC(iy,iz,it,iylp,Tee,offset,1) = sum( WYm(:,:,:,:,iylp,Tee,1) * actionTz(:,:,:,:) ) / volume
                           actionYpC(iy,iz,it,iylp,Tee,offset,2) = sum( WYp(:,:,:,:,iylp,Tee,2) * actionTz(:,:,:,:) ) / volume
                           actionYmC(iy,iz,it,iylp,Tee,offset,2) = sum( WYm(:,:,:,:,iylp,Tee,2) * actionTz(:,:,:,:) ) / volume

                           topChgYpC(iy,iz,it,iylp,Tee,offset,1) = sum( WYp(:,:,:,:,iylp,Tee,1) * topChgTz(:,:,:,:) ) / volume
                           topChgYmC(iy,iz,it,iylp,Tee,offset,1) = sum( WYm(:,:,:,:,iylp,Tee,1) * topChgTz(:,:,:,:) ) / volume
                           topChgYpC(iy,iz,it,iylp,Tee,offset,2) = sum( WYp(:,:,:,:,iylp,Tee,2) * topChgTz(:,:,:,:) ) / volume
                           topChgYmC(iy,iz,it,iylp,Tee,offset,2) = sum( WYm(:,:,:,:,iylp,Tee,2) * topChgTz(:,:,:,:) ) / volume

                        end do  ! iylp loop end
                     end do     ! it = 0, nt-1 loop end
                  end do        ! iz loop end
               end do           ! iy loop end
            end do              ! Tee loop end
         end do                 ! offset loop end
         call system_clock(time2)
         write(*,'(a,f15.8,a)') 'time spent in the offset loop is: ',(real(time2-time1)/real(count_rate)),' s'
         filename=output(1:len_trim(output))//'.Yp-xy'//fstr1(1:len_trim(fstr1))
         call writeYshape(filename,actionYpC(:,:,:,:,:,:,1),WYpavg(:,:,1),avgAction,offmax)

         filename=output(1:len_trim(output))//'.Yp-xy'//fstr2(1:len_trim(fstr2))
         call writeYshape(filename,TopChgYpC(:,:,:,:,:,:,1),WYpavg(:,:,1),avgTopChg,offmax)

         filename=output(1:len_trim(output))//'.Ym-xy'//fstr1(1:len_trim(fstr1))
         call writeYshape(filename,actionYmC(:,:,:,:,:,:,1),WYmavg(:,:,1),avgAction,offmax)

         filename=output(1:len_trim(output))//'.Ym-xy'//fstr2(1:len_trim(fstr2))
         call writeYshape(filename,TopChgYpC(:,:,:,:,:,:,1),WYmavg(:,:,1),avgTopChg,offmax)

         filename=output(1:len_trim(output))//'.Yp-xz'//fstr1(1:len_trim(fstr1))
         call writeYshape(filename,actionYpC(:,:,:,:,:,:,2),WYpavg(:,:,2),avgAction,offmax)

         filename=output(1:len_trim(output))//'.Yp-xz'//fstr2(1:len_trim(fstr2))
         call writeYshape(filename,TopChgYpC(:,:,:,:,:,:,2),WYpavg(:,:,2),avgTopChg,offmax)

         filename=output(1:len_trim(output))//'.Ym-xz'//fstr1(1:len_trim(fstr1))
         call writeYshape(filename,actionYmC(:,:,:,:,:,:,2),WYmavg(:,:,2),avgAction,offmax)

         filename=output(1:len_trim(output))//'.Ym-xz'//fstr2(1:len_trim(fstr2))
         call writeYshape(filename,TopChgYpC(:,:,:,:,:,:,2),WYmavg(:,:,2),avgTopChg,offmax)

 999     CALL SYSTEM_CLOCK(end_count)

         elapsed_time = REAL(end_count - start_count) / REAL(count_rate)
c
c     Excluding assignmentgs
c
         if (elapsed_time .ne. 0.0) then

            write(*,'(/,a,f20.5,a,i3.3)') 'The elapsed time is',elapsed_time,
     &           ' seconds in configuration:',cfg

         endif

      end do                    !end cfg loop

      end program VacuumRespY

      

