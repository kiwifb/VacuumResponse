c
c     Program to illustrate correlations between static quarks in a baryon and 
c     gluon action or topological charge
c     Deals with with L shaped Wilson loops. The loops do not match Y loops
c     but the program design is inspired by the TY program.
c     Adapted from VacuumRespY_plan.f by FB050324
c     This program study L shapes that make isocelis triangles (iqx==iqy)
c     8 shapes are created, 4 in each of the 2 possible plans.
c     The results are averaged on a single L shape (know in other program as L1-xy).
c
      program VacuumRespLY

      include'VRfiles/VRCommonUse.h'
      USE L_LYLOOPSRX
      USE L_WRITEYSHAPE
      USE L_PRODUCT

      IMPLICIT NONE

c     general variables
      integer,parameter                                          :: offmax=4
      integer,parameter                                          :: nsteps=4
      double precision, parameter                                :: alpha=0.7d0
c
      integer                                                    :: offset,smear3d
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
!HPF$ TEMPLATE WLloops(nx,ny,nz,nt,0:nYloop,4)
!HPF$ DISTRIBUTE WLloops(BLOCK,BLOCK,*,*,*,*)
      double precision,dimension(nx,ny,nz,nt)                    :: reconAction
      double precision,dimension(nx,ny,nz,nt)                    :: topQDensity
!HPF$ ALIGN (:,:,*,*) WITH WLloops(:,:,*,*,*,*) :: reconAction,topQDensity
      double precision,dimension(nx,ny,nz,nt)                    :: actionTx,actionTy,actionTz
!HPF$ ALIGN (:,:,*,*) WITH WLloops(:,:,*,*,*,*) :: actionTx,actionTy,actionTz
      double precision,dimension(nL(that),offmax)                :: avgAction, avgTopChg
!HPF$ DISTRIBUTE avgAction(BLOCK,*)
!HPF$ DISTRIBUTE avgTopChg(BLOCK,*)
      double precision,dimension(nx,ny,nz,nt)                    :: topChgTx,topChgTy,topChgTz
!HPF$ ALIGN (:,:,*,*) WITH WLloops(:,:,*,*,*) :: topChgTx,topChgTy,topChgTz
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
      character(len=132)                                         :: filename
c
c  Time averaged Quantities
c
      double precision,dimension(nx,ny,nz,nt,0:nYloop,4)         :: WLxy,WLxz ! T-shape wilson loops
!HPF$ ALIGN (:,:,*,*,*,*) WITH WLloops(:,:,*,*,*,*) :: WLxy,WLxz

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

      double precision,dimension(0:nYloop,nL(that))              :: WLavg
!HPF$ DISTRIBUTE WLavg(BLOCK,*)
      double precision,dimension(nx,ny,nz,nt,nc,nc)              :: uptr,upti
!HPF$ DISTRIBUTE uptr(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE upti(*,*,BLOCK,BLOCK,*,*)
      double precision,dimension(nx,ny,nz,nt)                    :: actionT
      double precision,dimension(nx,ny,nz,nt)                    :: topChgT
!HPF$ ALIGN (:,:,*,*) WITH WLloops(:,:,*,*,*,*) :: actionT,topChgT
c
c  Correlations
c
      double precision,dimension(0:ny-1,0:nz-1,0:nt-1,0:nYloop,nL(that),offmax) :: actionLC
      double precision,dimension(0:ny-1,0:nz-1,0:nt-1,0:nYloop,nL(that),offmax) :: topChgLC
!HPF$ DISTRIBUTE actionLC(*,*,*,BLOCK,BLOCK,*) 
!HPF$ DISTRIBUTE topChgLC(*,*,*,BLOCK,BLOCK,*) 
c
c  Counting indexes
c
      integer                                                    :: iylp,iydir,ish
      integer                                                    :: ict,icy,icz

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
      write(*,'(a)')'Please enter the amount of spatial smearing.'
      read (*,*) smear3d
      write(*,*) smear3d

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

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c'     initialize the epsilon index for the calls from ly_loop                    !' 
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

            write (output,'(a,a,i3.3,2(a,i2.2))')trim(prefix),trim(lastconfig),cfg,'.s',nsteps,'-s',smear3d

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
c  We will put the electric field in the reconAction and the
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
c     
      actionLC = 0.0d0
      topChgLC = 0.0d0

      call system_clock(time1)

      do Tee = 1, nL(that)

         call product(utr,uti,that,uptr,upti,Tee)

         call ly_loops_mirror(xhat,yhat,utr,uti,uptr,upti,Tee,WLxy)
         call ly_loops_mirror(xhat,zhat,utr,uti,uptr,upti,Tee,WLxz)

         do iylp = 1, nYloop
                
            WLavg(iylp,Tee) = (Sum( WLxy(:,:,:,:,iylp,:)) + 
     &                         Sum( WLxz(:,:,:,:,iylp,:))) / (8*volume)

         end do

         write(*,*) Tee
         write(*,*) WLavg(:,Tee)

         do offset = 1, offmax

            if ( Tee - 2*offset .lt. 0 ) then
               cycle
            end if
             
            do it = 1,nx

               actionT(it,:,:,:) = 0.0d0
               topChgT(it,:,:,:) = 0.0d0

               do is = it + offset , it + Tee - offset
                     
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

                     do iylp = 1, nYloop

                        actionLC( iy, iz, it,iylp,Tee,offset) = actionLC( iy, iz, it,iylp,Tee,offset) +
     &                           sum( WLxy(:,:,:,:,iylp,1) * actionTz(:,:,:,:) )
                        
                        topChgLC( iy, iz, it,iylp,Tee,offset) = topChgLC( iy, iz, it,iylp,Tee,offset) +
     &                           sum( WLxy(:,:,:,:,iylp,1) * topChgTz(:,:,:,:) )
                        
                        ict = modulo(nt - it , nt)

                        actionLC( iy, iz,ict,iylp,Tee,offset) = actionLC( iy, iz,ict,iylp,Tee,offset) +
     &                           sum( WLxy(:,:,:,:,iylp,4) * actionTz(:,:,:,:) )
                        
                        topChgLC( iy, iz,ict,iylp,Tee,offset) = topChgLC( iy, iz,ict,iylp,Tee,offset) +
     &                           sum( WLxy(:,:,:,:,iylp,4) * topChgTz(:,:,:,:) )
                        
                        icz = modulo(nz - iz , nz)

                        actionLC( iy,icz, it,iylp,Tee,offset) = actionLC( iy,icz, it,iylp,Tee,offset) +
     &                           sum( WLxy(:,:,:,:,iylp,2) * actionTz(:,:,:,:) )
                        
                        topChgLC( iy,icz, it,iylp,Tee,offset) = topChgLC( iy,icz, it,iylp,Tee,offset) +
     &                           sum( WLxy(:,:,:,:,iylp,2) * topChgTz(:,:,:,:) )
                        
                        actionLC( iy,icz,ict,iylp,Tee,offset) = actionLC( iy,icz,ict,iylp,Tee,offset) +
     &                           sum( WLxy(:,:,:,:,iylp,3) * actionTz(:,:,:,:) )
                        
                        topChgLC( iy,icz,ict,iylp,Tee,offset) = topChgLC( iy,icz,ict,iylp,Tee,offset) +
     &                           sum( WLxy(:,:,:,:,iylp,3) * topChgTz(:,:,:,:) )
                        
                        icz = modulo(     iy , nz)
                        icy = modulo(ny - iz , ny)

                        actionLC(icy,icz, it,iylp,Tee,offset) = actionLC(icy,icz, it,iylp,Tee,offset) +
     &                           sum( WLxz(:,:,:,:,iylp,1) * actionTz(:,:,:,:) )
                        
                        topChgLC(icy,icz, it,iylp,Tee,offset) = topChgLC(icy,icz, it,iylp,Tee,offset) +
     &                           sum( WLxz(:,:,:,:,iylp,1) * topChgTz(:,:,:,:) )

                        ict  = modulo(nt - it , nt)

                        actionLC(icy,icz,ict,iylp,Tee,offset) = actionLC(icy,icz,ict,iylp,Tee,offset) +
     &                           sum( WLxz(:,:,:,:,iylp,4) * actionTz(:,:,:,:) )
                        
                        topChgLC(icy,icz,ict,iylp,Tee,offset) = topChgLC(icy,icz,ict,iylp,Tee,offset) +
     &                           sum( WLxz(:,:,:,:,iylp,4) * topChgTz(:,:,:,:) )

                        icz = modulo(nz - iy , nz)
                        icy = modulo(     iz , ny)

                        actionLC(icy,icz, it,iylp,Tee,offset) = actionLC(icy,icz, it,iylp,Tee,offset) +
     &                           sum( WLxz(:,:,:,:,iylp,2) * actionTz(:,:,:,:) )
                        
                        topChgLC(icy,icz, it,iylp,Tee,offset) = topChgLC(icy,icz, it,iylp,Tee,offset) +
     &                           sum( WLxz(:,:,:,:,iylp,2) * topChgTz(:,:,:,:) )

                        actionLC(icy,icz,ict,iylp,Tee,offset) = actionLC(icy,icz,ict,iylp,Tee,offset) +
     &                           sum( WLxz(:,:,:,:,iylp,3) * actionTz(:,:,:,:) )
                        
                        topChgLC(icy,icz, it,iylp,Tee,offset) = topChgLC(icy,icz, it,iylp,Tee,offset) +
     &                           sum( WLxz(:,:,:,:,iylp,3) * topChgTz(:,:,:,:) )


                     end do     ! iylp loop end
                  end do        ! it = 0, nt-1 loop end
               end do           ! iz loop end
            end do              ! iy loop end
         end do                 ! offset loop end
      end do                    ! Tee loop end

      call system_clock(time2)
      write(*,'(a,f15.8,a)') 'time spent in the offset loop is: ',(real(time2-time1)/real(count_rate)),' s'

      actionLC(:,:,:,1:nYloop,:,:) = actionLC(:,:,:,1:nYloop,:,:) / (8 * volume)
      topChgLC(:,:,:,1:nYloop,:,:) = topChgLC(:,:,:,1:nYloop,:,:) / (8 * volume)
      
c      do iylp = 1,nYloop
         
c         actionLC(:,:,:,iylp,:,:) = cshift(cshift(
c     &        actionLC(:,:,:,iylp,:,:),dim=3,shift=iylp/2),dim=2,shift=iylp/2)
c         actionLC(:,:,:,iylp,:,:) = cshift(cshift(
c     &        actionLC(:,:,:,iylp,:,:),dim=3,shift=iylp/2),dim=2,shift=iylp/2)

c         topChgLC(:,:,:,iylp,:,:) = cshift(cshift(
c     &        topChgLC(:,:,:,iylp,:,:),dim=3,shift=iylp/2),dim=2,shift=iylp/2)
c         topChgLC(:,:,:,iylp,:,:) = cshift(cshift(
c     &        topChgLC(:,:,:,iylp,:,:),dim=3,shift=iylp/2),dim=2,shift=iylp/2) 

c      end do
      
      write(filename,'(a,a,a)') trim(output),'.Lavg8',trim(fstr1)
      call writeYshape(filename,actionLC(:,:,:,:,:,:),WLavg(:,:),avgAction,offmax)

      write(filename,'(a,a,a)') trim(output),'.Lavg8',trim(fstr2)
      call writeYshape(filename,TopChgLC(:,:,:,:,:,:),WLavg(:,:),avgTopChg,offmax)

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

      end program VacuumRespLY
      
