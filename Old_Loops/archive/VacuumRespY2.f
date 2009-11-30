c
c     Program to illustrate correlations between static quarks in a baryon and 
c     gluon action or topological charge
c     Deals with with Y shaped Wilson loops.
c     Adapted from VacuumRespLT.f by F. Bissey, Jan. 2004
c
      program VacuumRespY2

      USE baryonParam
      USE TADPOLEIMP
      USE READLINKS
      USE ape_smear
      USE ape_smear3D
      USE fMuNu
      USE TopQandReconAction
      USE EleMagField
      USE YLOOPSSELECT

      IMPLICIT NONE

c     general variables
      integer,parameter                                          :: offmax=4
      integer,parameter                                          :: nsteps=4
      integer,parameter                                          :: nYloop=4
      double precision, parameter                                :: alpha=0.7d0
c
      integer                                                    :: offset
      double precision                                           :: beta     
      logical                                                    :: smear_links
      double precision,dimension(nx,ny,nz,nt,mu,nc,nc)           :: ur,ui,utr,uti
!HPF$ DISTRIBUTE ur(*,*,BLOCK,BLOCK,*,*,*)
!HPF$ DISTRIBUTE ui(*,*,BLOCK,BLOCK,*,*,*)
!HPF$ DISTRIBUTE utr(*,*,BLOCK,BLOCK,*,*,*)
!HPF$ DISTRIBUTE uti(*,*,BLOCK,BLOCK,*,*,*)


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
      integer                                                    :: nxold,nyold,nzold,ntold
      double precision                                           :: betaold
      double precision,dimension(nx,ny,nz,nt)                    :: action
!HPF$ DISTRIBUTE action(*,*,BLOCK,BLOCK)
      double precision,dimension(nx,ny,nz,nt,nc,nc,nf)           :: Fr,Fi        !Fmunu
!HPF$ DISTRIBUTE Fr(*,*,BLOCK,BLOCK,*,*,*)
!HPF$ DISTRIBUTE Fi(*,*,BLOCK,BLOCK,*,*,*)
      double precision,dimension(nx,ny,nz,nt)                    :: reconAction
!HPF$ DISTRIBUTE reconAction(*,*,BLOCK,BLOCK)
      double precision,dimension(nx,ny,nz,nt)                    :: topQDensity
!HPF$ DISTRIBUTE topQDensity(*,*,BLOCK,BLOCK)
      double precision                                           :: Q=0.0d0 
c
c  Time averaged Quantities
c
      double precision,dimension(nx,ny,nz,nt,0:nYloop,1:nL(that))   :: WY ! Y-shape wilson loops
!HPF$ DISTRIBUTE WY(*,*,BLOCK,BLOCK,*,*)

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Nomenclature for the Y shape Wilson loop
c     The index running from 0 to nY correspond to the various
c     position of the three quarks in the Y-loop.
c     Here are the sets of coordinates for each value of the index:
c
c     0 : all 3 quarks at (0,0)
c     1 :
c     2 :
c     3 :
c     4 :
c     5 :
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      double precision,dimension(0:nYloop,1:nL(that))            :: WYavg
!HPF$ DISTRIBUTE WYavg(*,*)

      double precision,dimension(nx,ny,nz,nt)                    :: actionT
!HPF$ DISTRIBUTE actionT(*,*,BLOCK,BLOCK)
      double precision,dimension(nx,ny,nz,nt)                    :: actionTx,actionTy,actionTz
!HPF$ DISTRIBUTE actionTx(*,*,BLOCK,BLOCK)
!HPF$ DISTRIBUTE actionTy(*,*,BLOCK,BLOCK)
!HPF$ DISTRIBUTE actionTz(*,*,BLOCK,BLOCK)
      double precision,dimension(nL(that),offmax)                :: avgAction, avgTopChg
!HPF$ DISTRIBUTE avgAction(*,*)
!HPF$ DISTRIBUTE avgTopChg(*,*)
      double precision                                           :: maxAction, maxTopChg
      double precision                                           :: minAction, minTopChg
      double precision,dimension(nx,ny,nz,nt)                    :: topChgT
!HPF$ DISTRIBUTE topChgT(*,*,BLOCK,BLOCK)
      double precision,dimension(nx,ny,nz,nt)                    :: topChgTx,topChgTy,topChgTz
!HPF$ DISTRIBUTE topChgTx(*,*,BLOCK,BLOCK)
!HPF$ DISTRIBUTE topChgTy(*,*,BLOCK,BLOCK)
!HPF$ DISTRIBUTE topChgTz(*,*,BLOCK,BLOCK)
c
c  Correlations
c
      double precision,dimension(0:ny-1,0:nz-1,0:nt-1,0:nYloop,1:nL(that),offmax) :: actionYC
!HPF$ DISTRIBUTE actionYC(*,BLOCK,BLOCK,*,*,*)
      double precision,dimension(0:ny-1,0:nz-1,0:nt-1,0:nYloop,1:nL(that),offmax) :: topChgYC
!HPF$ DISTRIBUTE topChgYC(*,BLOCK,BLOCK,*,*,*)
c
c  Counting indexes
c
      integer                                                    :: iy,iz,it,ic,Tee,is,iiy
c
c  Timer support
c
      integer                                                    :: start_count, end_count, count_rate
      integer                                                    :: time1, time2
      real                                                       :: elapsed_time

c
c  Begin execution
c
      CALL SYSTEM_CLOCK(start_count, count_rate)
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

      basedir = basedir(1:len_trim(basedir))//lastconfig(1:len_trim(lastconfig))

      call ReadLinks(basedir,ur,ui,nfig,beta,nxold,nyold,nzold,ntold,lastPlaq,plaqbarAvg,uzero)

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
     &          ( 8.0d0*pi**2/(6.0d0/beta) )

        reconAction = reconAction * beta * (mu*(mu-1)/2.0d0) / ( 8.0d0*pi**2/(6.0d0/beta) )

        write(*,'(/,2(a,f15.8,/),/)') 
     &                'Reconstructed S/S_0 = ', sons1,
     &                'Topological Charge  = ', Q

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      else

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  We'll put the electric field in the reconAction and the
c                magnetic field in the topQdensity
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

         call EleMagField(Fr,Fi,uzero,reconAction,topQdensity)

         write(*,'(/,2(a,e15.8,/),/)') 
     &                'Sum of Electric Field = ', SUM(reconAction), 
     &                'Sum of Magnetic Field = ', SUM(topQdensity)

      end if

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  Get the Wilson Loops
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      utr = ur
      uti = ui

      do istep = 1,10
            
         call ape_smear3D(utr,uti,alpha,that)

      end do

      if ( smear_links ) then

         write (string,'(a,i2.0)') '.s', nsteps
         lastConfig = prefix(1:len_trim(prefix))//lastConfig(1:len_trim(lastConfig))//string

      else

         lastConfig = prefix(1:len_trim(prefix))//lastConfig(1:len_trim(lastConfig))//'.s00'

      endif

      if (Correlate == 1) then

         open(11,file=lastConfig(1:len_trim(lastConfig))//'.Y.action.correl.unf',status='unknown',form='unformatted')

         open(12,file=lastConfig(1:len_trim(lastConfig))//'.Y.topChg.correl.unf',status='unknown',form='unformatted')

      else

         open(11,file=lastConfig(1:len_trim(lastConfig))//'.Y.ele.correl.unf',status='unknown',form='unformatted')

         open(12,file=lastConfig(1:len_trim(lastConfig))//'.Y.mag.correl.unf',status='unknown',form='unformatted')

      end if

      call y_loops(ur,ui,WY)

      do iiy = 0, nYloop

         do it = 2, nL(that)
               
            WYavg(iiy,it) = Sum( WY(:,:,:,:,iiy,it) ) / volume

         end do                 ! it = 2, nL(that)  loop end
      end do

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

c     write(*,'(a,i2,a,i2)') 'Step ',iy,' of ',ny-1
               
               actionTx(:,:,:,:) = cshift(actionT(:,:,:,:), dim = 2, shift = iy)
               topChgTx(:,:,:,:) = cshift(topChgT(:,:,:,:), dim = 2, shift = iy)

               do iz = 0, nz-1

                  actionTy(:,:,:,:) = cshift(actionTx(:,:,:,:), dim = 3, shift = iz)
                  topChgTy(:,:,:,:) = cshift(topChgTx(:,:,:,:), dim = 3, shift = iz)

                  do it = 0, nt-1
                     
                     actionTz(:,:,:,:) = cshift(actionTy(:,:,:,:), dim = 4, shift = it)
                     topChgTz(:,:,:,:) = cshift(topChgTy(:,:,:,:), dim = 4, shift = it)

                     do iiy = 0, nYloop

                        actionYC(iy,iz,it,iiy,Tee,offset) = sum( WY(:,:,:,:,iiy,Tee) * actionTz(:,:,:,:) ) / volume

                        topChgYC(iy,iz,it,iiy,Tee,offset) = sum( WY(:,:,:,:,iiy,Tee) * topChgTz(:,:,:,:) ) / volume

                     end do     ! iiy loop end
                  end do        ! it = 0, nt-1 loop end
               end do           ! iz loop end
            end do              ! iy loop end
         end do                 ! Tee loop end
      end do                    ! offset loop end
      call system_clock(time2)
      write(*,'(a,f15.8,a)') 'time spent in the offset loop is: ',(real(time2-time1)/real(count_rate)),' s'

      write(11) actionYC(:,:,:,:,:,:)
      write(11) WYavg(:,:), avgAction(:,:)

      write(12) topChgYC(:,:,:,:,:,:)
      write(12) WYavg(:,:), avgTopChg(:,:)
            
      close(11)
      close(12)

 999  CALL SYSTEM_CLOCK(end_count)

      elapsed_time = REAL(end_count - start_count) / REAL(count_rate)
c
c  Excluding assignmentgs
c
      if (elapsed_time .ne. 0.0) then

         write(*,'(/,a,f20.5,a)') 'The elapsed time is',elapsed_time,
     &   ' seconds for program VacuumResp.'

      endif

      end program VacuumRespY2

      

