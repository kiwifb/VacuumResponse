c
c  A program to average the results from VacuumRespY-plan
c  Adapted from AverageY
c
c  v.1.01 AK040622
c  v.1.02 FB040624 converted DOS-MAC to unix line end
c  v.1.03 FB040624 corrected and improved conversion various check
c  v.1.04 FB040625 the bugsquashing sessions
c  v.1.1  FB041201 Accept both Y and TY shapes and writes correct filenames
c  v.1.2  FB050209 Modified to accept (T)Yavg4 and extract WL at the same time
c         or perform potential analisys from potential.f 
c  v.1.3  FB050223 Modified to allow "averaging by part" 
c  v.2.0  FB050322 Make the potential analyses a Jack-knife one.
c  v.2.0a FB050323 Corrected _the_one_bug_ I found from the previous day modifications
c                  Correction may trigger cleanup later on. 
c  v.2.1  FB050415 Various corrections since last version bump.
c                  Split in a legacy version (2.0a) and the new version
c                  dropping old stuff. This program can now deal with
c                  Y, TY and LY loops. It also can spit out some DX output
c                  but it is limited to avg# kind of data.
c  v.2.2  FB050520 Finalized centering of LY shapes and DX output.
c  v.2.3  FB050808 Corrected enormous bug in the potential code. It was
c                  harmless if you had 200 configurations.
c  v.2.4  FB051028 Added an output option for measuring the radius of the node.
c  v.3.0  FB060314 Adding support for computing errors on each points of the lattice-
c                  Cleaning variable names as well.
c  v.3.0a FB060317 Cleanning and adjusting all the code modifications.
c  v.3.0.9 FB060727 Adding support for the new "DiQuark" configuration.
c  v.3.1.0 FB060727 Finished support for DiQuark configurations. Rebuilt the
c                    change lost by importing the doublehelix tree from the JK
c                    fork. Some cleaning and restructuring.
c  v.3.2.0 FB061009 Added QQbar support.
c  v.3.2.1 FB061123 Properly center topological charge + correcting various blunders.
c                    that didn't affect the "action" results.
c
      program AverageYall
c
      USE L_baryonParam
      USE L_WRITEAVG
c
      IMPLICIT NONE
c
c     for s16t32
      integer,parameter                                       :: maxsize=12
c     for s12t24:
c     integer,parameter                                       :: maxsize=9
      integer                                                 :: nYloop
      integer,parameter                                       :: offmax=4
      integer                                                 :: start, finish, icon, Correlate
      integer                                                 :: iiy,offset,Tee,icon2
      integer                                                 :: ix,iy,iz,sx,sy
      integer                                                 :: shpe,zero
      integer                                                 :: pot,dx,nodex,nodey,xerr
      character(len=80)                                       :: configName, suffixName, reportName
      character(len=132)                                      :: thisConfig,dxname
      character(len=12)                                       :: fstr1,fstr2,wstr1,wstr2
      character(len=30)                                       :: dialog
c
      double precision                                        :: maxAction,  maxTopChg
      double precision                                        :: minAction,  minTopChg
      double precision                                        :: maxActionPref,  maxTopChgPref
      double precision                                        :: minActionPref,  minTopChgPref
c
      double precision,dimension(nL(that),offmax)             :: avgAction, avgTopChg
!HPF$ DISTRIBUTE avgAction(BLOCK,BLOCK)
!HPF$ DISTRIBUTE avgTopChg(BLOCK,BLOCK)
      double precision,dimension(nL(that),offmax)             :: avgActionA, avgTopChgA
!HPF$ DISTRIBUTE avgActionA(BLOCK,BLOCK)
!HPF$ DISTRIBUTE avgTopChgA(BLOCK,BLOCK)
c
      double precision,dimension(0:ny-1,0:nz-1,0:nt-1,0:maxsize,nL(that),offmax)   :: actionC
!HPF$ DISTRIBUTE actionC(*,BLOCK,BLOCK,*,*,*)
      double precision,dimension(0:ny-1,0:nz-1,0:nt-1,0:maxsize,nL(that),offmax)   :: topChgC
!HPF$ DISTRIBUTE topChgC(*,BLOCK,BLOCK,*,*,*)
      double precision,dimension(0:ny-1,0:nz-1,0:nt-1,0:maxsize,nL(that),offmax)   :: actionA
!HPF$ DISTRIBUTE actionA(*,BLOCK,BLOCK,*,*,*)
      double precision,dimension(0:ny-1,0:nz-1,0:nt-1,0:maxsize,nL(that),offmax)   :: topChgA
!HPF$ DISTRIBUTE topChgA(*,BLOCK,BLOCK,*,*,*)
      double precision,dimension(0:maxsize,nL(that))          :: Wavg
!HPF$ DISTRIBUTE Wavg(BLOCK,BLOCK)
      double precision,dimension(0:maxsize,nL(that))          :: WavgA
!HPF$ DISTRIBUTE WavgA(BLOCK,BLOCK)
      double precision,dimension(0:ny-1,0:nz-1,0:nt-1,0:maxsize,nL(that),offmax)   :: CAction
      double precision,dimension(0:ny-1,0:nz-1,0:nt-1,0:maxsize,nL(that),offmax)   :: CTopChg
      double precision,dimension(0:ny-1,0:nz-1,0:nt-1,0:maxsize,nL(that),offmax) :: ErrCAction
      double precision,dimension(0:maxsize,nL(that),200)      :: WavgG,Vi
      double precision,dimension(0:maxsize,nL(that)-1)        :: V0,Vb,Vberr
      double precision,dimension(maxsize)                     :: rs, dqq
      double precision,dimension(0:nt-1,0:nz-1,0:ny-1)        :: dxac,dxtc
      integer,dimension(maxsize)                              :: shx
c
c  Begin execution
c
c  Basic input
c
      write(*,*)
      write(*,'(a)')'What would you like to correlate?'
      write(*,'(a)')'      1: Action and topological charge.'
      write(*,'(a)')'      2: Electric and Magnetic fields.'
      read (*,*) Correlate
c
      write(*,'(a)') 'What is the base file name?'
      read (*,'(a80)') configName
c
      write(*,'(a)') 'Please provide a range of configuration numbers. e.g. 1 100'
      read (*,*) start, finish
c
c  The suffix is of the form ".s 4" for example.
c
      write(*,'(a)') 'What is the file name suffix?'
      read (*,'(a80)') suffixName
c
      write(*,'(a)') 'Are we dealing with:'
      write(*,'(a)') '1) Y shape'
      write(*,'(a)') '2) TY shape'
      write(*,'(a)') '3) LY shape'
      write(*,'(a)') '4) DQ2 shape'
      write(*,'(a)') '5) DQ4 shape'
      write(*,'(a)') '6) QQbar shape'
      write(*,'(a)') 'results?'
      read (*,*) shpe
c
c  --------------
c  average distance to the Fermat point (rs)
c  average quark separation (dqq)
c  --------------
c
      select case (shpe)
      case(:2)

         nYloop = 9

         rs(1) = 0.15301
         rs(2) = 0.26503
         rs(3) = 0.30603
         rs(4) = 0.41804
         rs(5) = 0.57106
         rs(6) = 0.72407
         rs(7) = 0.87708
         rs(8) = 0.9891
         rs(9) = 1.14211

         dqq(1) = 0.25568
         dqq(2) = 0.39383
         dqq(3) = 0.51136
         dqq(4) = 0.64907
         dqq(5) = 0.90455
         dqq(6) = 1.16012
         dqq(7) = 1.41573
         dqq(8) = 1.55359
         dqq(9) = 1.80911

      case (3)

         nYloop = 11

         rs(1) = 0.07921
         rs(2) = 0.15841
         rs(3) = 0.23762
         rs(4) = 0.31682
         rs(5) = 0.39603
         rs(6) = 0.47524
         rs(7) = 0.55444
         rs(8) = 0.63365
         rs(9) = 0.71285
         rs(10) = 0.79206
         rs(11) = 0.87127

         dqq(1) = 0.13998
         dqq(2) = 0.27997
         dqq(3) = 0.41995
         dqq(4) = 0.55993
         dqq(5) = 0.69991
         dqq(6) = 0.83990
         dqq(7) = 0.97988
         dqq(8) = 1.11986
         dqq(9) = 1.25984
         dqq(10) = 1.39983
         dqq(11) = 1.53981

      case (4)

         nYloop = 12

         do iiy = 1, nYloop

            rs(iiy)  = (iiy+sqrt(3.d0))*0.123/3.d0
            dqq(iiy) = 2.d0*(1.d0+sqrt(iiy**2+1.d0))*0.123/3.d0

         end do

      case (5)

         nYloop = 12

         do iiy = 1, nYloop

            rs(iiy)  = (iiy+2.d0*sqrt(3.d0))*0.123/3.d0
            dqq(iiy) = 2.d0*(2.d0+sqrt(iiy**2+4.d0))*0.123/3.d0

         end do

      case (6)

         nYloop = 12

            do iiy = 1, nYloop

               rs(iiy) = iiy*0.123d0
               dqq(iiy) = rs(iiy)

            end do

      end select

      write(*,*) 'nYloop=',nYloop
c
c  ------------
c  data for the baseline record
c  ------------
c
      shx(1) = 1
      shx(2) = 1
      shx(3) = 1
      shx(4) = 2
      shx(5) = 3
      shx(6) = 4
      shx(7) = 4
      shx(8) = 4
      shx(9) = 5
c
c  -----------
c  Options
c  -----------
c
      write(*,'(a)') 'Do you want a potential study with that (yes=1 no=0)?'
      read (*,*) pot
c
      write(*,'(a)') 'Do you want an error to be computed for each lattice point (yes=1 no=0)?'
      write(*,'(a)') '(Only for the correlated action for now.)'
      read (*,*) xerr
c
      write(*,'(a)') 'Do you want a openDX output (yes=1 no=0)?'
      read (*,*) dx
c
      write(*,'(a)') 'Do you want a baseline output to measure the size of the node  (yes=1 no=0)?'
      read (*,*) nodex
c
      write(*,'(a)') 'Do you want a baseline output to measure the size of the flux tube  (yes=1 no=0)?'
      read (*,*) nodey
c
      if( (pot+xerr+dx+nodex+nodey) == 0 ) then 
         write(*,'(a)') 'are you averaging by part (yes=1 no=0)?'
         read (*,*) zero
      else
         zero = 0
      end if
c
      write(*,'(a)') 'Please enter a report file name.'
      read (*,'(a80)') reportName
c
c  Initialize accumulators
c
      ErrCAction   = 0.d0
      actionA      = 0.0d0
      topChgA      = 0.0d0
      WavgA        = 0.0d0
      avgActionA   = 0.0d0
      avgTopChgA   = 0.0d0
      WavgG        = 0.0d0
c
c Initialiaze file strings and dialogue strings
c
      if ( Correlate == 1) then
c
         fstr1='.action'
         fstr2='.topChg'
         wstr1='Action'
         wstr2='TopChg'
         dialog='action and topological '
c
      else
c
         fstr1='.ele'
         fstr2='.mag'
         wstr1='Electric'
         wstr2='Magnetic'
         dialog='electric and magnetic field '
c
      endif
c
      select case (shpe)
      case (1)
         write(suffixName,fmt='(2a)') trim(suffixName),'.Yavg4'
         if(zero == 0) then
            write(reportName,fmt='(2a)') trim(reportName),'.Y'
         end if

      case (2)
         write(suffixName,fmt='(2a)') trim(suffixName),'.TYavg4'
         if(zero == 0) then
            write(reportName,fmt='(2a)') trim(reportName),'.TY'
         end if
      case (3)
         write(suffixName,fmt='(2a)') trim(suffixName),'.Lavg8'
         if(zero == 0) then
            write(reportName,fmt='(2a)') trim(reportName),'.L'
         end if
      case (4)
         write(suffixName,fmt='(2a)') trim(suffixName),'.DQ2avg4'
         if(zero == 0) then
            write(reportName,fmt='(2a)') trim(reportName),'.DQ2'
         end if
      case (5)
         write(suffixName,fmt='(2a)') trim(suffixName),'.DQ4avg4'
         if(zero == 0) then
            write(reportName,fmt='(2a)') trim(reportName),'.DQ4'
         end if
      case (6)
         write(suffixName,fmt='(2a)') trim(suffixName),'.QQ'
         if(zero == 0) then
            write(reportName,fmt='(2a)') trim(reportName),'.QQ'
         end if
      end select
c
      do icon = start, finish
c
c  Create file name
c
         write(thisConfig,fmt='(a,i3.3,a)') trim(configName),icon,trim(suffixName)
c
         write(*,*) trim(thisConfig)
c
c  -------------------------------------
c  read in results for one configuration
c  -------------------------------------
c
         write(*,'(3a)') 'Reading ',dialog,'correlations.'
c
         open (11,file=trim(thisConfig)//trim(fstr1)//'.correl.unf',status='old',form='unformatted')
         open (12,file=trim(thisConfig)//trim(fstr2)//'.correl.unf',status='old',form='unformatted')
c
c  ----------
c  accumulate
c  ----------
c
c   avg4/avg8
c
         read(11) actionC(:,:,:,0:nYloop,:,:)
         read(11) Wavg(0:nYloop,:), avgAction(:,:)
         read(12) topChgC(:,:,:,0:nYloop,:,:)
         read(12) Wavg(0:nYloop,:), avgTopChg(:,:)
c
         ActionA    = ActionA    + ActionC
         topChgA    = topChgA    + topChgC
         avgActionA = avgActionA + avgAction
         avgTopChgA = avgTopChgA + avgTopChg
         WavgA      = WavgA      + Wavg

         close(11)
         close(12)
c
c----------------------------------------
c potential computations
c----------------------------------------
c
         if((pot == 1).or.(xerr == 1)) then

            do icon2 = start, finish

               if(icon2 /= icon) then
                  WavgG(1:nYloop,:,icon2) = WavgG(1:nYloop,:,icon2) + Wavg(1:nYloop,:)
               end if

            end do

         end if
c
      end do                    ! end icon loop
c
c  -------
c  average
c  -------
c
      actionA       = actionA    / (finish - start + 1 )
      topChgA       = topChgA    / (finish - start + 1 )
      WavgA         = WavgA      / (finish - start + 1 )
      avgActionA    = avgActionA / (finish - start + 1 )
      avgTopChgA    = avgTopChgA / (finish - start + 1 )


      if((pot == 1).or.(xerr == 1)) then
         WavgG      = WavgG      / (finish - start)
      end if

      if(zero==0) then ! We are not doing an averaging by part.
c
c  --------------
c  Computing the correlated action
c  --------------
c
         do Tee = 2,nL(that)
            do offset = 1,offmax
c
               if ( Tee - 2*offset .lt. 0 ) then
                  cycle
               endif
c
               do iiy = 1, nYloop
c
                  CAction(:,:,:,iiy,Tee,offset) =
     &                 actionA(:,:,:,iiy,Tee,offset) / (WavgA(iiy,Tee) * avgActionA(Tee,offset))
c
                  CTopChg(:,:,:,iiy,Tee,offset) =
     &                 topChgA(:,:,:,iiy,Tee,offset) / (WavgA(iiy,Tee) * avgTopChgA(Tee,offset))

               end do
            end do
         end do
c
c  --------------
c  Errors on the correlated action
c  --------------
c
         if(xerr == 1) then

            do icon = start, finish

               write(thisConfig,fmt='(a,i3.3,a)') trim(configName),icon,trim(suffixName)
               open (11,file=trim(thisConfig)//trim(fstr1)//'.correl.unf',status='old',form='unformatted')

               read(11) actionC(:,:,:,0:nYloop,:,:)
               read(11) Wavg(0:nYloop,:), avgAction(:,:)

               close(11)

               do Tee = 2,nL(that)
                  do offset = 1,offmax
c
                     if ( Tee - 2*offset .lt. 0 ) then
                        cycle
                     end if
c
                     do iiy = 1, nYloop
c
                        ErrCAction(:,:,:,iiy,Tee,offset) = ErrCAction(:,:,:,iiy,Tee,offset) +
     &                      (CAction(:,:,:,iiy,Tee,offset) -
     &                         actionC(:,:,:,iiy,Tee,offset)/(Wavg(iiy,Tee)*avgAction(Tee,offset)))**2

                     end do ! iiy loop
                  end do    ! offset loop
               end do       ! Tee loop
            end do          ! icon loop

            ErrCAction = sqrt(ErrCAction / ((finish - start + 1)*(finish - start)))

         end if
c
c  --------------------------------
c  move L, DQ  and QQbar shapes into the center
c  --------------------------------
c
c  L shape
c
         if(shpe == 3) then

            do iiy = 1,nYloop

               actionA(:,:,:,iiy,:,:) = cshift(
     &              cshift(actionA(:,:,:,iiy,:,:), shift=iiy/2,dim=3),
     &              shift=iiy/2,dim=2)

               CAction(:,:,:,iiy,:,:) = cshift(
     &              cshift(CAction(:,:,:,iiy,:,:), shift=iiy/2,dim=3),
     &              shift=iiy/2,dim=2)

               if(xerr == 1) then

                  ErrCAction(:,:,:,iiy,:,:) = cshift(
     &                 cshift(ErrCAction(:,:,:,iiy,:,:), shift=iiy/2,dim=3),
     &                 shift=iiy/2,dim=2)

               end if

               topChgA(:,:,:,iiy,:,:) = cshift(
     &              cshift(topChgA(:,:,:,iiy,:,:), shift=iiy/2,dim=3),
     &              shift=iiy/2,dim=2)

               CTopChg(:,:,:,iiy,:,:) = cshift(
     &              cshift(CTopChg(:,:,:,iiy,:,:), shift=iiy/2,dim=3),
     &              shift=iiy/2,dim=2)

            end do

         end if
c
c  DQ shapes && QQbar shape
c
         if((shpe.eq.4).or.(shpe.eq.5).or.(shpe.eq.6)) then

            do iiy = 1,nYloop

               actionA(:,:,:,iiy,:,:) =
     &              cshift(actionA(:,:,:,iiy,:,:), shift=ny/2,dim=2)

               CAction(:,:,:,iiy,:,:) =
     &              cshift(CAction(:,:,:,iiy,:,:), shift=ny/2,dim=2)

               if(xerr == 1) then

                  ErrCAction(:,:,:,iiy,:,:) =
     &                 cshift(ErrCAction(:,:,:,iiy,:,:), shift=ny/2,dim=2)

               end if

               topChgA(:,:,:,iiy,:,:) =
     &              cshift(topChgA(:,:,:,iiy,:,:), shift=ny/2,dim=2)

               CTopChg(:,:,:,iiy,:,:) =
     &              cshift(CTopChg(:,:,:,iiy,:,:), shift=ny/2,dim=2)

            end do

         end if
c
c ------------------------------------------
c baseline record before centering.
c ------------------------------------------
c
c ------------------------------------------
c baseline record: node to quark in the long direction.
c ------------------------------------------
c

         if(nodex == 1) then

            do iiy = 1, nYloop

               if(xerr == 1) then
                  write(thisconfig,fmt='(3a,i2.2,a)') 'BaseX/',trim(reportName),'-node.',iiy,'-err.dat'
               else
                  write(thisconfig,fmt='(3a,i2.2,a)') 'BaseX/',trim(reportName),'-node.',iiy,'.dat'
               end if

               open (11,file=trim(thisconfig),status='unknown',form='formatted')

               do ix =-nt/2, nt/2-1

                  if(shpe==2) then
                     sx = modulo(nt+ix-shx(iiy),nt)
                  else
                     sx = modulo(nt+ix,nt)
                  end if

                  if(xerr == 1) then
                     write(11,*) ix,CAction(0,0,sx,iiy,2:nL(that),1),ErrCAction(0,0,sx,iiy,2:nL(that),1)
                  else
                     write(11,*) ix,CAction(0,0,sx,iiy,2:nL(that),1)
                  end if

               end do

               close(11)

            end do

         end if
c
c ------------------------------------------
c baseline record: flux tube diameter measurements.
c ------------------------------------------
c
      if(nodey==1) then

         do iiy = 1, nYloop

            do ix= -nt/2, nt/2-1

               if(shpe==2) then
                  sx = modulo(nt+ix-shx(iiy),nt)
               else
                  sx = modulo(nt+ix,nt)
               end if

               if(xerr == 1) then
                  write(thisconfig,fmt='(3a,i2.2,a,i3.2,a)') 'BaseY/',trim(reportName),'-flux.',iiy,'x',ix,'-err.dat'
               else
                  write(thisconfig,fmt='(3a,i2.2,a,i3.2,a)') 'BaseY/',trim(reportName),'-flux.',iiy,'x',ix,'.dat'
               end if

               open (11,file=trim(thisconfig),status='unknown',form='formatted')
c
               do iy =-ny/2, ny/2-1

                  sy = modulo(ny+iy,ny)

                  if(xerr == 1) then
                     write(11,*) iy,CAction(0,sy,sx,iiy,2:nL(that),1),ErrCAction(0,sy,sx,iiy,2:nL(that),1)
                  else
                     write(11,*) iy,CAction(0,sy,sx,iiy,2:nL(that),1)
                  end if

               end do

               close(11)

            end do

         end do

      end if
c
c  --------------------------------
c  move origin to center of lattice
c  --------------------------------
c
         actionA(:,:,:,:,:,:) = cshift(
     &        cshift(
     &        cshift(actionA(:,:,:,:,:,:), shift=-nt/2, dim=3),
     &        shift=-nz/2, dim=2),
     &        shift=-ny/2, dim=1)

         topChgA(:,:,:,:,:,:) = cshift(
     &        cshift(
     &        cshift(topChgA(:,:,:,:,:,:), shift=-nt/2, dim=3),
     &        shift=-nz/2, dim=2),
     &        shift=-ny/2, dim=1)

         CAction(:,:,:,:,:,:) = cshift(
     &        cshift(
     &        cshift(CAction(:,:,:,:,:,:), shift=-nt/2, dim=3),
     &        shift=-nz/2, dim=2),
     &        shift=-ny/2, dim=1)

         CTopChg(:,:,:,:,:,:) = cshift(
     &        cshift(
     &        cshift(CTopChg(:,:,:,:,:,:), shift=-nt/2, dim=3),
     &        shift=-nz/2, dim=2),
     &        shift=-ny/2, dim=1)

         if(xerr == 1) then

            ErrCAction(:,:,:,:,:,:) = cshift(
     &           cshift(
     &           cshift(CAction(:,:,:,:,:,:), shift=-nt/2, dim=3),
     &           shift=-nz/2, dim=2),
     &           shift=-ny/2, dim=1)


         end if
c
c----------------------------------
c potential computations
c----------------------------------
c
         if(pot == 1) then

            do icon = start, finish

               do Tee = 1, nL(that) - 1

                  do iiy = 1, nYloop

                     Vi(iiy,Tee,icon) = log(WavgG(iiy,Tee,icon)/WavgG(iiy,Tee+1,icon))

                  end do
               end do

            end do

            do Tee = 1, nL(that) - 1
               do iiy = 1, nYloop

                  V0(iiy,Tee) = log(WavgA(iiy,Tee)/WavgA(iiy,Tee+1))

               end do
            end do

            do iiy = 1, nYloop
               do Tee = 1, nL(that) - 1

                  Vb(iiy,Tee) = sum(Vi(iiy,Tee,start:finish)) /(finish - start +1)
                  Vberr(iiy,Tee) = sqrt(
     &                 ((sum((Vi(iiy,Tee,start:finish)-Vb(iiy,Tee))**2))/
     &                 (finish - start +1))*(finish - start))

               end do
            end do
c
c----------------------------------
c record potential results
c----------------------------------
c
            select case (shpe)
            case (1)
               write(suffixName,fmt='(a)') '.wilson-Y.dat'
            case (2)
               write(suffixName,fmt='(a)') '.wilson-T.dat'
            case (3)
               write(suffixName,fmt='(a)') '.wilson-L.dat'
            case (4)
               write(suffixName,fmt='(a)') '.wilson-DQ2.dat'
            case (5)
               write(suffixName,fmt='(a)') '.wilson-DQ4.dat'
            case (6)
               write(suffixName,fmt='(a)') '.wilson-QQ.dat'
            end select

            do iiy = 1, nYloop

               write(thisconfig,fmt='(3a,i2.2,a)') 'Pot/',trim(reportName),'-sh',iiy,trim(suffixName) 

               open (11,file=trim(thisconfig),status='unknown',form='formatted')
c
               do Tee = 1, nL(that) - 1

                  if(WavgA(iiy,Tee+1)<= 0.d0) then 
                     exit
                  end if

                  write (11,*) V0(iiy,Tee),Vberr(iiy,Tee),Vb(iiy,Tee)

               end do

               close(11)

            end do

            write(thisconfig,fmt='(4a)') 'Pot/',trim(reportName),'-pot',trim(suffixName) 

            open (11,file=trim(thisconfig),status='unknown',form='formatted')

            do iiy = 1, nYloop

               write(11,*) rs(iiy),dqq(iiy),V0(iiy,1),Vberr(iiy,1),Vb(iiy,1)

            end do

            close(11)

            write(thisconfig,fmt='(4a)') 'Pot/',trim(reportName),'-pot-T2',trim(suffixName) 

            open (11,file=trim(thisconfig),status='unknown',form='formatted')

            do iiy = 1, nYloop

               if(WavgA(iiy,3)<0.d0) then 
                  cycle
               end if

               write(11,*) rs(iiy),dqq(iiy),V0(iiy,2),Vberr(iiy,2),Vb(iiy,2)

            end do

            close(11)

            write(thisconfig,fmt='(4a)') 'Pot/',trim(reportName),'-pot-T3',trim(suffixName) 

            open (11,file=trim(thisconfig),status='unknown',form='formatted')

            do iiy = 1, nYloop

               if(WavgA(iiy,4)<0.d0) then
                  cycle
               end if

                write(11,*) rs(iiy),dqq(iiy),V0(iiy,3),Vberr(iiy,3),Vb(iiy,3)

            end do

            close(11)

         end if
c
c  --------------
c  record results
c  --------------
c
c  --------------
c  record the averaged data set
c  --------------
c
         open (9,file=trim(reportName)//'.report',status='unknown')
c
         open (11,file='Averages/'//trim(reportName)//trim(fstr1)//'.correl.unf',status='unknown',form='unformatted')
         open (12,file='Averages/'//trim(reportName)//trim(fstr2)//'.correl.unf',status='unknown',form='unformatted')
c
         do Tee = 2,nL(that)
            do offset = 1,offmax
c
               if ( Tee - 2*offset .lt. 0 ) then
                  cycle
               endif
c
c     This part of the program has been changed to be more like the writing in the AverageT
c     program and easily allow the use of unf_writeavg. Previous version of the program didn't
c     have the loop over iiy and the data for 12t24 are still in old fashion writing at this
c     stage (18 Oct. 2004).
c
               do iiy = 1,nYloop

                  call unf_writeavg(11,actionA(:,:,:,iiy,Tee,offset),WavgA(iiy,Tee),avgActionA(Tee,offset))
                  call unf_writeavg(12,topChgA(:,:,:,iiy,Tee,offset),WavgA(iiy,Tee),avgtopChgA(Tee,offset))

               end do
c
            enddo                  ! end offset loop
         enddo                     ! end Tee loop
c
         close(11)
         close(12)
c
c  In AVS format unscaled
c
         open (11,file='Averages/'//trim(reportName)//trim(fstr1)//'.correl',status='unknown')
         open (12,file='Averages/'//trim(reportName)//trim(fstr2)//'.correl',status='unknown')
         if(xerr == 1) then
            open (13,file='Averages/'//trim(reportName)//trim(fstr1)//'.err-correl',status='unknown')
         end if
c
         do Tee = 2,nL(that)
            do offset = 1,offmax
c
               if ( Tee - 2*offset .lt. 0 ) then
                  cycle
               endif
c
               do iiy = 1, nYloop
c
                  write(11,'(a)')    '*********************'
                  write(11,'(a,3i3)') 'Quark separation scheme ', iiy, Tee,offset
                  write(11,'(a)')    '*********************'
                  call fmt_writeavg(11,CAction(:,:,:,iiy,Tee,offset))
c
                  if(xerr == 1) then
                     write(13,'(a)')    '*********************'
                     write(13,'(a,3i3)') 'Quark separation scheme ', iiy, Tee,offset
                     write(13,'(a)')    '*********************'
                     call fmt_writeavg(13,ErrCAction(:,:,:,iiy,Tee,offset))
                  end if

                  write(12,'(a)')    '*********************'
                  write(12,'(a,3i3)') 'Quark separation scheme ', iiy, Tee,offset
                  write(12,'(a)')    '*********************'
                  call fmt_writeavg(12,CTopChg(:,:,:,iiy,Tee,offset))
c
               enddo               ! end iiy loop
c
            enddo                  ! end offset loop
         enddo                     ! end Tee loop
c
         close (11)
         close (12)
         if(xerr == 1) then
            close(13)
         end if
c
c  Normalized to 0 ... 255
c
         open (11,file='Averages/'//trim(reportName)//trim(fstr1)//'.correl.255',status='unknown')
         open (12,file='Averages/'//trim(reportName)//trim(fstr2)//'.correl.255',status='unknown')
         if(xerr == 1) then
            open (13,file='Averages/'//trim(reportName)//trim(fstr1)//'.err-correl.255',status='unknown')
         end if
c
c  Trap preferred values
c
         do Tee = 2,nL(that)
            do offset = 1,offmax
c
               if ( Tee - 2*offset .lt. 0 ) then
                  cycle
               endif
c
               maxActionPref = 1.0082
               minActionPref = 0.9118
               maxTopChgPref = 1.0174
               minTopChgPref = 0.8478
c
               do iiy = 1, nYloop
c
c  report max min values of correlations
c     for action and topological charge
c
                  write(9,'(/,a,3i3,/,4(3a,f8.4,/))') 
     &                 'At separation scheme, ', iiy,Tee,offset,
     &                 'Max',wstr1,'Correlation = ', maxval(CAction(:,:,:,iiy,Tee,offset)),
     &                 'Min',wstr1,'Correlation = ', minval(CAction(:,:,:,iiy,Tee,offset)),
     &                 'Max',wstr2,'Correlation = ', maxval(CTopChg(:,:,:,iiy,Tee,offset)),
     &                 'Min',wstr2,'Correlation = ', minval(CTopChg(:,:,:,iiy,Tee,offset))
c
                  maxAction = max(maxval(CAction(:,:,:,iiy,Tee,offset)),maxActionPref)
                  minAction = min(minval(CAction(:,:,:,iiy,Tee,offset)),minActionPref)
                  CAction(:,:,:,iiy,Tee,offset) = (CAction(:,:,:,iiy,Tee,offset)-minAction) * 255.0d0 / 
     &                 (maxAction - minAction)
c
                  if(xerr == 1) then
                     ErrCAction(:,:,:,iiy,Tee,offset) = ErrCAction(:,:,:,iiy,Tee,offset) * 255.0d0 / 
     &                    (maxAction - minAction)
                  end if

                  maxTopChg = max(maxval(CTopChg(:,:,:,iiy,Tee,offset)),maxTopChgPref)
                  minTopChg = min(minval(CTopChg(:,:,:,iiy,Tee,offset)),minTopChgPref)
                  CTopChg(:,:,:,iiy,Tee,offset) = (CTopChg(:,:,:,iiy,Tee,offset)-minTopChg) * 255.0d0 / 
     &                 (maxTopChg-minTopChg)
c
                  write(11,'(a)')    '***********************'
                  write(11,'(a,3i3)') 'Quark separation scheme ', iiy, Tee, offset
                  write(11,'(a)')    '***********************'
                  call fmt_writeavg(11,CAction(:,:,:,iiy,Tee,offset))
c
                  if(xerr == 1) then
                     write(13,'(a)')    '***********************'
                     write(13,'(a,3i3)') 'Quark separation scheme ', iiy, Tee, offset
                     write(13,'(a)')    '***********************'
                     call fmt_writeavg(13,ErrCAction(:,:,:,iiy,Tee,offset))
                  end if

                  write(12,'(a)')    '***********************'
                  write(12,'(a,3i3)') 'Quark separation scheme', iiy, Tee,offset
                  write(12,'(a)')    '***********************'
                  call fmt_writeavg(12,CTopChg(:,:,:,iiy,Tee,offset))

                  if(dx == 1) then

                     write(dxname,'(4a,3(i2.2))') 'DX-out/',trim(reportName),trim(fstr1),'.correl.DX-',iiy,Tee,offset

                     open (14,file=dxname,status='unknown')

                     write(dxname,'(4a,3(i2.2))') 'DX-out/',trim(reportName),trim(fstr2),'.correl.DX-',iiy,Tee,offset

                     open (15,file=dxname,status='unknown')

                     do ix = 0, nt-1
                        do iy = 0, nz-1
                           do iz = 0, ny-1

                              dxac(ix,iy,iz) = CAction(iz,iy,ix,iiy,Tee,offset)
                              dxtc(ix,iy,iz) = CTopChg(iz,iy,ix,iiy,Tee,offset)

                           end do
                        end do
                     end do

                     call fmt_writeavg(14,dxac(:,:,:))
                     call fmt_writeavg(15,dxtc(:,:,:))

                     close(14)
                     close(15)

                  end if
c
               enddo               ! end iiy loop
c
            enddo                  ! end offset loop
         enddo                     ! end Tee loop
c
         close (11)
         close (12)
         if(xerr == 1) then
            close(13)
         end if
c
         close (9)
c
      else ! Averaging by parts - recording the partial average.
c
         open (11,file=trim(reportName)//trim(fstr1)//'.correl.unf',status='unknown',form='unformatted')
         open (12,file=trim(reportName)//trim(fstr2)//'.correl.unf',status='unknown',form='unformatted')
c
         write(11) ActionA(:,:,:,0:nYloop,:,:)
         write(11) WavgA(0:nYloop,:), avgActionA(0:nYloop,:)
         write(12) topChgA(:,:,:,0:nYloop,:,:)
         write(12) WavgA(0:nYloop,:), avgtopChgA(0:nYloop,:)
c
         close(11)
         close(12)
c
      endif
c
      end program AverageYall
