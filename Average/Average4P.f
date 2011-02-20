c
c  A program to average the results for 4 quarks wilson loops
c  Branched from UniversalAverage on 20100513
c  The main difference is that some loops are 2 dimensional.
c
c  v.0.1  FB20100513 Initial work on this incarnation of the program.
c
      program Average4P
c
      USE L_baryonParam
      USE L_WRITEAVG
      USE A_TRAPEZE
c
      IMPLICIT NONE
c
c     for s16t32
      integer,parameter                                       :: mxsize=12
      integer,parameter                                       :: mysize=6
c     for s12t24:
c     integer,parameter                                       :: maxsize=9
      integer                                                 :: nlx,nly
      integer,parameter                                       :: offmax=4
      integer                                                 :: start, finish, icon, Correlate
      integer                                                 :: ilx,ily,offset,Tee,icon2
      integer                                                 :: ix,iy,iz,sx,sy
      integer                                                 :: shpe,zero
      integer                                                 :: pot,dx,nodex,nodey,xerr,vtk
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
      double precision,dimension(0:ny-1,0:nz-1,0:nt-1,1:mxsize,1:mysize,nL(that),offmax)   :: actionC
!HPF$ DISTRIBUTE actionC(*,BLOCK,BLOCK,*,*,*,*)
      double precision,dimension(0:ny-1,0:nz-1,0:nt-1,1:mxsize,1:mysize,nL(that),offmax)   :: topChgC
!HPF$ DISTRIBUTE topChgC(*,BLOCK,BLOCK,*,*,*,*)
      double precision,dimension(0:ny-1,0:nz-1,0:nt-1,1:mxsize,1:mysize,nL(that),offmax)   :: actionA
!HPF$ DISTRIBUTE actionA(*,BLOCK,BLOCK,*,*,*,*)
      double precision,dimension(0:ny-1,0:nz-1,0:nt-1,1:mxsize,1:mysize,nL(that),offmax)   :: topChgA
!HPF$ DISTRIBUTE topChgA(*,BLOCK,BLOCK,*,*,*,*)
      double precision,dimension(1:mxsize,1:mysize,nL(that))          :: Wavg
!HPF$ DISTRIBUTE Wavg(*,BLOCK,BLOCK)
      double precision,dimension(1:mxsize,1:mysize,nL(that))          :: WavgA
!HPF$ DISTRIBUTE WavgA(*,BLOCK,BLOCK)
      double precision,dimension(0:ny-1,0:nz-1,0:nt-1,1:mxsize,1:mysize,nL(that),offmax)   :: CAction
      double precision,dimension(0:ny-1,0:nz-1,0:nt-1,1:mxsize,1:mysize,nL(that),offmax)   :: CTopChg
      double precision,dimension(0:ny-1,0:nz-1,0:nt-1,1:mxsize,1:mysize,nL(that),offmax) :: ErrCAction
      double precision,dimension(1:mxsize,1:mysize,nL(that),200)      :: WavgG,Vi
      double precision,dimension(1:mxsize,1:mysize,nL(that)-1)        :: V0,Vb,Vberr
      double precision,dimension(1:mxsize,1:mysize)                     :: rs, dqq
      double precision                                        :: lsize
      double precision,dimension(0:nt-1,0:nz-1,0:ny-1)        :: dxac,dxtc
      integer,dimension(1:mxsize,1:mysize)                              :: shx
      integer,dimension(nL(that),1:mxsize,1:mysize)                     :: exept
      integer,dimension(nL(that),1:mxsize,1:mysize,200)                 :: flag
      integer                                                 :: expuls
      double precision,dimension(1:mxsize,1:mysize,nL(that),200)      :: Cexpuls
      double precision,dimension(1:mxsize,1:mysize,nL(that))          :: VExp,VExpErr
      double precision                                        :: tmpExpuls
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
      write(*,'(a)') '1) DDQ-PL shape'
      write(*,'(a)') '2) DDQ-OP shape'
      write(*,'(a)') '3) 4P shape'
      write(*,'(a)') 'results?'
      read (*,*) shpe
c
c  Ask for the lattice spacing
c
      write(*,'(a)') 'What is the lattice spacing in fermi? (1 if you do not know)'
      read (*,*) lsize
c
c  --------------
c  average distance to the Fermat point (rs)
c  average quark separation (dqq)
c  --------------
c
      select case (shpe)
      case(:2)

         nlx= 12
         nly= 6

         do ilx = 1, nlx
            do ily = 1, nly

               rs(ilx,ily)=
              dqq(ilx,ily)= (4.d0*ily+2.d0*ilx)*lsize
            end do
         end do

      case(3)

         nlx= 6

         do ilx = 1, nlx

           rs(ilx,1)=

      case (5)

         nloop = 11

         do iiy = 1, nloop

            rs(iiy)  = (1.d0+sqrt(3.d0))*iiy*lsize/(3.d0*sqrt(2.d0))
            dqq(iiy) = (2.d0+sqrt(2.d0))*iiy*lsize/3.d0

         end do

      case (6)

         nloop = 12

         do iiy = 1, nloop

            rs(iiy)  = (iiy+sqrt(3.d0))*lsize/3.d0
            dqq(iiy) = 2.d0*(1.d0+sqrt(iiy**2+1.d0))*lsize/3.d0

         end do

      case (7)

         nloop = 12

         do iiy = 1, nloop

            rs(iiy)  = (iiy+2.d0*sqrt(3.d0))*lsize/3.d0
            dqq(iiy) = 2.d0*(2.d0+sqrt(iiy**2+4.d0))*lsize/3.d0

         end do

      case (8)

         nloop = 12

            do iiy = 1, nloop

               rs(iiy) = iiy*lsize
               dqq(iiy) = rs(iiy)

            end do

      end select

      write(*,*) 'nloop=',nloop
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
      write(*,'(a)') 'Do you want a vtk output (yes=1 no=0)?'
      read (*,*) vtk
c
      write(*,'(a)') 'Do you want a baseline output to measure the size of the node  (yes=1 no=0)?'
      read (*,*) nodex
c
      write(*,'(a)') 'Do you want a baseline output to measure the size of the flux tube  (yes=1 no=0)?'
      read (*,*) nodey
c
      write(*,'(a)') 'Do you want a report on the expulsion (yes=1 no=0)?'
      read (*,*) expuls
c
      if( (pot+xerr+dx+nodex+nodey+expuls+vtk) == 0 ) then
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
      ErrCAction   = 0.0d0
      actionA      = 0.0d0
      topChgA      = 0.0d0
      WavgA        = 0.0d0
      avgActionA   = 0.0d0
      avgTopChgA   = 0.0d0
      WavgG        = 0.0d0

      Cexpuls      = 0.0d0
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
         write(suffixName,fmt='(2a)') trim(suffixName),'.VYavg4'
         if(zero == 0) then
            write(reportName,fmt='(2a)') trim(reportName),'.VY'
         end if

      case (4)
         write(suffixName,fmt='(2a)') trim(suffixName),'.VY1avg4'
         if(zero == 0) then
            write(reportName,fmt='(2a)') trim(reportName),'.VY1'
         end if
      case (5)
         write(suffixName,fmt='(2a)') trim(suffixName),'.Lavg8'
         if(zero == 0) then
            write(reportName,fmt='(2a)') trim(reportName),'.L'
         end if
      case (6)
         write(suffixName,fmt='(2a)') trim(suffixName),'.DQ2avg4'
         if(zero == 0) then
            write(reportName,fmt='(2a)') trim(reportName),'.DQ2'
         end if
      case (7)
         write(suffixName,fmt='(2a)') trim(suffixName),'.DQ4avg4'
         if(zero == 0) then
            write(reportName,fmt='(2a)') trim(reportName),'.DQ4'
         end if
      case (8)
         write(suffixName,fmt='(2a)') trim(suffixName),'.QQ'
         if(zero == 0) then
            write(reportName,fmt='(2a)') trim(reportName),'.QQ'
         end if
      case (9)
         write(suffixName,fmt='(2a)') trim(suffixName),'.DLTavg4'
         if(zero == 0) then
            write(reportName,fmt='(2a)') trim(reportName),'.DLTavg4'
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
         read(11) actionC(:,:,:,0:nloop,:,:)
         read(11) Wavg(0:nloop,:), avgAction(:,:)
         read(12) topChgC(:,:,:,0:nloop,:,:)
         read(12) Wavg(0:nloop,:), avgTopChg(:,:)
c
         ActionA    = ActionA    + ActionC
         topChgA    = topChgA    + topChgC
         avgActionA = avgActionA + avgAction
         avgTopChgA = avgTopChgA + avgTopChg
         WavgA      = WavgA      + Wavg

         close(11)
         close(12)
c
c  --------------
c  Computing the correlated action for the expulsion calculation
c  --------------
c
         if(expuls == 1) then

            do Tee = 2,nL(that)
c  Only offset = 1 for now
               offset = 1
c
               do iiy = 1, nloop
c  1-C for the expulsion amount
                  CAction(:,:,:,iiy,Tee,offset) = 1.0d0 -
     &                 actionC(:,:,:,iiy,Tee,offset) / (Wavg(iiy,Tee) * avgAction(Tee,offset))
c
               tmpExpuls = trapeze3D(CAction(:,:,:,iiy,Tee,offset),ny,nz,nt)/(ny*nz*nt)
               do icon2 = start, finish

                  if(icon2 /= icon) then
                     Cexpuls(iiy,Tee,icon2) = Cexpuls(iiy,Tee,icon2) + tmpExpuls
                  end if

               end do

               end do
            end do

         endif
c
c----------------------------------------
c potential computations
c----------------------------------------
c
         if(pot == 1) then

            do icon2 = start, finish

               if(icon2 /= icon) then
                  WavgG(1:nloop,:,icon2) = WavgG(1:nloop,:,icon2) + Wavg(1:nloop,:)
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


      if(pot == 1) then
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
               do iiy = 1, nloop
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

               read(11) actionC(:,:,:,0:nloop,:,:)
               read(11) Wavg(0:nloop,:), avgAction(:,:)

               close(11)

               do Tee = 2,nL(that)
                  do offset = 1,offmax
c
                     if ( Tee - 2*offset .lt. 0 ) then
                        cycle
                     end if
c
                     do iiy = 1, nloop
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
         if(shpe == 5) then

            do iiy = 1,nloop

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
         if((shpe.eq.6).or.(shpe.eq.7).or.(shpe.eq.8)) then

            do iiy = 1,nloop

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

            do iiy = 1, nloop

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

         do iiy = 1, nloop

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

            exept = 0
            flag  = 0

            do icon = start, finish

               do Tee = 1, nL(that) - 1

                  do iiy = 1, nloop

                     if((WavgG(iiy,Tee,icon)/WavgG(iiy,Tee+1,icon))>0.d0) then
                        Vi(iiy,Tee,icon) = log(WavgG(iiy,Tee,icon)/WavgG(iiy,Tee+1,icon))
                     else
                        Vi(iiy,Tee,icon) = 0.d0
                        exept(Tee,iiy) = exept(Tee,iiy) + 1
                        flag(Tee,iiy,icon) = 1
                     end if

                  end do
               end do

            end do

            do Tee = 1, nL(that) - 1
               do iiy = 1, nloop

                  if((WavgA(iiy,Tee)/WavgA(iiy,Tee+1))>0.d0) then
                     V0(iiy,Tee) = log(WavgA(iiy,Tee)/WavgA(iiy,Tee+1))
                  else
                     V0(iiy,Tee) = 0.d0
                  end if

               end do
            end do

            do iiy = 1, nloop
               do Tee = 1, nL(that) - 1

                  Vb(iiy,Tee) = sum(Vi(iiy,Tee,start:finish)) /(finish - start +1 - exept(Tee,iiy))

                  do icon = start, finish

                     if(flag(Tee,iiy,icon) == 1) then
                        Vi(iiy,Tee,icon) = Vb(iiy,Tee)
                     end if

                  end do

                  if((finish - start +1 - exept(Tee,iiy)) == 0) then
                     Vberr(iiy,Tee) = 0.d0
                  else
                     Vberr(iiy,Tee) = sqrt(
     &                    ((sum((Vi(iiy,Tee,start:finish)-Vb(iiy,Tee))**2))/
     &                    (finish - start +1 - exept(Tee,iiy)))*(finish - start - exept(Tee,iiy)))
                  end if

               end do
            end do
c
c----------------------------------
c record potential results
c----------------------------------
c
            write(suffixName,fmt='(a)') '.wilson.dat'

            do iiy = 1, nloop

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

            do iiy = 1, nloop

               write(11,*) rs(iiy),dqq(iiy),V0(iiy,:),Vberr(iiy,:)

            end do

            close(11)

         end if
c
c----------------------------------
c Expulsion computations
c----------------------------------
c
         if(expuls == 1) then

            do iiy = 1, nloop
               do Tee = 1, nL(that) - 1

                  VExp(iiy,Tee) = sum(Cexpuls(iiy,Tee,start:finish)) /(finish - start +1)

                  VExpErr(iiy,Tee) = sqrt(
     &                    ((sum((Cexpuls(iiy,Tee,start:finish)-VExp(iiy,Tee))**2))/
     &                    (finish - start +1))*(finish - start))

               end do
            end do
c
c----------------------------------
c record Expulsion results
c----------------------------------
c
            write(suffixName,fmt='(a)') '.expulsion.dat'

            do Tee = 1, nL(that) - 1

               write(thisconfig,fmt='(3a,i2.2,a)') 'Pot/',trim(reportName),'-Tee',Tee,trim(suffixName)

               open (11,file=trim(thisconfig),status='unknown',form='formatted')
c
               do iiy = 1, nloop
                  write (11,*) iiy,VExp(iiy,Tee),VExpErr(iiy,Tee)
               end do

               close(11)

            end do

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
               do iiy = 1,nloop

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
        if(vtk == 1) then
            open (14,file='VTK/'//trim(reportName)//trim(fstr1)//'.vtk',status='unknown')
            open (15,file='VTK/'//trim(reportName)//trim(fstr2)//'.vtk',status='unknown')

            write(14,'(a)')     '# vtk DataFile Version 3.0'
            write(14,'(a)')     trim(reportName)//'\n'
            write(14,'(a)')     'ASCII'
            write(14,'(a)')     'DATASET STRUCTURED_POINTS'
            write(14,'(a,3i2)') 'DIMENSIONS ',nt,nz,ny
            write(14,'(a)')     'ORIGIN 0 0 0'
            write(14,'(a)')     'SPACING 1 1 1'
            write(14,'(a,i4)')  'POINT_DATA',nt*nz*ny

            write(15,'(a)')     '# vtk DataFile Version 3.0'
            write(15,'(a)')     trim(reportName)//'\n'
            write(15,'(a)')     'ASCII'
            write(15,'(a)')     'DATASET STRUCTURED_POINTS'
            write(15,'(a,3i2)') 'DIMENSIONS ',nt,nz,ny
            write(15,'(a)')     'ORIGIN 0 0 0'
            write(15,'(a)')     'SPACING 1 1 1'
            write(15,'(a,i4)')  'POINT_DATA',nt*nz*ny
        end if
c
         do Tee = 2,nL(that)
            do offset = 1,offmax
c
               if ( Tee - 2*offset .lt. 0 ) then
                  cycle
               endif
c
               do iiy = 1, nloop
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
                  if(vtk == 1) then
                    write(14,'(3(a,i),a)')'SCALARS action_shape_',iiy,'-Tee_',Tee,'-offset_',offset,' float'
                    write(14,'(a)')     'LOOKUP_TABLE default'

                    write(15,'(3(a,i),a)')'SCALARS TopologicalCharge_shape_',iiy,'-Tee_',Tee,'-offset_',offset,' float'
                    write(15,'(a)')     'LOOKUP_TABLE default'

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
         if(vtk == 1) then
            close(14)
            close(15)
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
               do iiy = 1, nloop
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
         write(11) ActionA(:,:,:,0:nloop,:,:)
         write(11) WavgA(0:nloop,:), avgActionA(0:nloop,:)
         write(12) topChgA(:,:,:,0:nloop,:,:)
         write(12) WavgA(0:nloop,:), avgtopChgA(0:nloop,:)
c
         close(11)
         close(12)
c
      endif
c
      end program UniversalAverage
