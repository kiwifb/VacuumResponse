c
c  A program to average the results from VacuumRespLT
c  Adapted from the original Average for ConstituentQuarks
c  by F. Bissey Jan. 2004
c
      program AverageLT
c
      USE L_baryonParam
c
      IMPLICIT NONE
c
      include'Yfiles/yLoopSize.h'
      integer,parameter                                       :: offmax=4
      integer                                                 :: start, finish, icon, Correlate, flag
      integer                                                 :: iiy,offset,Tee
      character(len=3)                                        :: cfg
      character(len=80)                                       :: configName, suffixName, reportName
      character(len=132)                                      :: thisConfig
      character(len=12)                                       :: fstr1,fstr2,wstr1,wstr2
      character(len=30)                                       :: dialog
c
      double precision                                        :: maxAction,  maxTopChg
      double precision                                        :: minAction,  minTopChg
      double precision                                        :: maxActionPref,  maxTopChgPref
      double precision                                        :: minActionPref,  minTopChgPref

      double precision,dimension(nL(that),offmax)             :: avgAction, avgTopChg
!HPF$ DISTRIBUTE avgAction(*,*)
!HPF$ DISTRIBUTE avgTopChg(*,*)
      double precision,dimension(nL(that),offmax)             :: avgActionA, avgTopChgA
!HPF$ DISTRIBUTE avgActionA(*,*)
!HPF$ DISTRIBUTE avgTopChgA(*,*)

      double precision,dimension(0:ny-1,0:nz-1,0:nt-1,0:nYloop,nL(that),offmax)   :: actionC
!HPF$ DISTRIBUTE actionC(*,BLOCK,BLOCK,*,*,*)
      double precision,dimension(0:ny-1,0:nz-1,0:nt-1,0:nYloop,nL(that),offmax)   :: topChgC
!HPF$ DISTRIBUTE topChgC(*,BLOCK,BLOCK,*,*,*)
      double precision,dimension(0:ny-1,0:nz-1,0:nt-1,0:nYloop,nL(that),offmax)   :: actionA
!HPF$ DISTRIBUTE actionA(*,BLOCK,BLOCK,*,*,*,*)
      double precision,dimension(0:ny-1,0:nz-1,0:nt-1,0:nYloop,nL(that),offmax)   :: topChgA
!HPF$ DISTRIBUTE topChgA(*,BLOCK,BLOCK,*,*,*,*)
      double precision,dimension(0:nYloop,nL(that))           :: Wavg
!HPF$ DISTRIBUTE Wavg(*,*)
      double precision,dimension(0:nYloop,nL(that))           :: WavgA
!HPF$ DISTRIBUTE WavgA(*,*)

c
c  Correlations
c
      write(*,*)
      write(*,'(a)')'What would you like to correlate?'
      write(*,'(a)')'      1: Action and topological charge.'
      write(*,'(a)')'      2: Electric and Magnetic fields.'
      read (*,*) Correlate

      write(*,'(a)') 'What is the base file name?'
      read (*,'(a80)') configName
c
      write(*,'(a)') 'Please provide a range of configuration numbers. e.g. 0 100'
      read (*,*) start, finish
c
c  The suffix is of the form ".s 4" for example.
c
      write(*,'(a)') 'What is the file name suffix?'
      read (*,'(a80)') suffixName
c
      write(*,'(a)') 'Please enter a report file name.'
      read (*,'(a80)') reportName
c
c  Initialize accumulators
c
      actionA      = 0.0d0
      topChgA      = 0.0d0
      WavgA        = 0.0d0
      avgActionA   = 0.0d0
      avgTopChgA   = 0.0d0
c
c Initialiaze file strings and dialogue strings
c
      if ( Correlate == 1) then

         fstr1='.action'
         fstr2='.topChg'
         wstr1='Action'
         wstr2='TopChg'
         dialog='action and topological '

      else

         fstr1='.ele'
         fstr2='.mag'
         wstr1='Electric'
         wstr2='Magnetic'
         dialog='electric and magnetic field '

      endif

      do icon = start, finish
c
c  Create file name
c
         write(thisConfig,fmt='(a,i3.3,a)') trim(configName),icon,trim(suffixName)
         
         write(*,*) trim(thisConfig)
c
c  Read in results for one configuration
c
         write(*,'(/,3a,/)') 'Reading ',dialog,'correlations.'
         open (11,file=trim(thisConfig)//'.Y'//trim(fstr1)//'.correl.unf',status='old',form='unformatted')
         open (12,file=trim(thisConfig)//'.Y'//trim(fstr2)//'.correl.unf',status='old',form='unformatted')

         read(11) actionC(:,:,:,:,:,:)
         read(11) Wavg(:,:), avgAction(:,:)
               
         read(12) topChgC(:,:,:,:,:,:)
         read(12) Wavg(:,:), avgTopChg(:,:)
c
c     Accumulate
c
         ActionA    = ActionA    + ActionC
         topChgA    = topChgA    + topChgC
         avgActionA = avgActionA + avgAction
         avgTopChgA = avgTopChgA + avgTopChg
         WavgA      = WavgA      + Wavg
         
         
         close(11)
         close(12)
           
      end do                    ! end icon loop
c
c  Average
c
      actionA      = actionA      / (finish - start + 1)
      topChgA      = topChgA      / (finish - start + 1)
      WavgA        = WavgA        / (finish - start + 1)
      avgActionA   = avgActionA   / (finish - start + 1)
      avgTopChgA   = avgTopChgA   / (finish - start + 1)

      actionA(:,:,:,:,:,:) = cshift(
     &     cshift(
     &     cshift(actionA(:,:,:,:,:,:), shift=-nt/2, dim=3),
     &     shift=-nz/2, dim=2),
     &     shift=-ny/2, dim=1)

      topChgA(:,:,:,:,:,:) = cshift(
     &     cshift(
     &     cshift(topChgA(:,:,:,:,:,:), shift=-nt/2, dim=3),
     &     shift=-nz/2, dim=2),
     &     shift=-ny/2, dim=1)


      open (9,file=trim(reportName)//'.report',status='unknown')

      open (11,file=trim(reportName)//'.Y'//trim(fstr1)//'.correl.unf',status='unknown',form='unformatted')
      open (12,file=trim(reportName)//'.Y'//trim(fstr2)//'.correl.unf',status='unknown',form='unformatted')

      do Tee = 2,nL(that)
         do offset = 1,offmax
 
            if ( Tee - 2*offset .lt. 0 ) then
               cycle
            endif

            write(11) actionA(:,:,:,:,Tee,offset)
            write(11) WavgA(:,Tee), avgActionA(Tee,offset)
            write(12) topChgA(:,:,:,:,Tee,offset)
            write(12) WavgA(:,Tee), avgTopChgA(Tee,offset)

         enddo                  ! end offset loop
      enddo                     ! end Tee loop

      close(11)
      close(12)
c
c  In AVS format unscaled
c
      open (11,file=trim(reportName)//'.Y'//trim(fstr1)//'.correl',status='unknown')
      open (12,file=trim(reportName)//'.Y'//trim(fstr2)//'.correl',status='unknown')
c
      do Tee = 2,nL(that)
         do offset = 1,offmax
            
            if ( Tee - 2*offset .lt. 0 ) then
               cycle
            endif
            
            do iiy = 0, nYloop

               write(11,'(a)')    '*********************'
               write(11,'(a,3i3)') 'Y Loop ', iiy, Tee,offset
               write(11,'(a)')    '*********************'
               write(11,'(f15.8)') actionA(:,:,:,iiy,Tee,offset) / (WavgA(iiy,Tee) * avgActionA(Tee,offset))

               write(12,'(a)')    '*********************'
               write(12,'(a,3i3)') 'Y Loop ', iiy, Tee,offset
               write(12,'(a)')    '*********************'
               write(12,'(f15.8)') topChgA(:,:,:,iiy,Tee,offset) / (WavgA(iiy,Tee) * avgTopChgA(Tee,offset))
            
            enddo               ! end iiy loop
            
         enddo                  ! end offset loop
      enddo                     ! end Tee loop

      close (11)
      close (12)
c
c  Normalized to 0 ... 255
c
      open (11,file=trim(reportName)//'.Y'//trim(fstr1)//'.correl.255',status='unknown')
      open (12,file=trim(reportName)//'.Y'//trim(fstr2)//'.correl.255',status='unknown')
c
c  Trap preferred values
c
      do Tee = 2,nL(that)
         do offset = 1,offmax

            if ( Tee - 2*offset .lt. 0 ) then
               cycle
            endif
             
c     maxActionPref = maxval(actionA(:,:,:,3,Tee,offset))
c     minActionPref = minval(actionA(:,:,:,3,Tee,offset))
c     maxTopChgPref = maxval(topChgA(:,:,:,3,Tee,offset))
c     minTopChgPref = minval(topChgA(:,:,:,3,Tee,offset))

            maxActionPref = 1.0082
            minActionPref = 0.9118
            maxTopChgPref = 1.0174
            minTopChgPref = 0.8478

            do iiy = 0, nYloop

               actionA(:,:,:,iiy,Tee,offset) = 
     &              actionA(:,:,:,iiy,Tee,offset) / (WavgA(iiy,Tee) * avgActionA(Tee,offset))
                   
               topChgA(:,:,:,iiy,Tee,offset) = 
     &              topChgA(:,:,:,iiy,Tee,offset) / (WavgA(iiy,Tee) * avgTopChgA(Tee,offset))
c
c  report max min values of correlations
c
               write(9,'(/,a,3i3,/,4(3a,f8.4,/))') 
     &              'At separation scheme, ', iiy,Tee,offset,
     &              'Max',wstr1,'Correlation = ', maxval(actionA(:,:,:,iiy,Tee,offset)),
     &              'Min',wstr1,'Correlation = ', minval(actionA(:,:,:,iiy,Tee,offset)),
     &              'Max',wstr2,'Correlation = ', maxval(topChgA(:,:,:,iiy,Tee,offset)),
     &              'Min',wstr2,'Correlation = ', minval(topChgA(:,:,:,iiy,Tee,offset))

               maxAction = max(maxval(actionA(:,:,:,iiy,Tee,offset)),maxActionPref)
               minAction = min(minval(actionA(:,:,:,iiy,Tee,offset)),minActionPref)
               actionA(:,:,:,iiy,Tee,offset) = (actionA(:,:,:,iiy,Tee,offset)-minAction) * 255.0d0 / 
     &              (maxAction - minAction)

               maxTopChg = max(maxval(topChgA(:,:,:,iiy,Tee,offset)),maxTopChgPref)
               minTopChg = min(minval(topChgA(:,:,:,iiy,Tee,offset)),minTopChgPref)
               topChgA(:,:,:,iiy,Tee,offset) = (topChgA(:,:,:,iiy,Tee,offset)-minTopChg) * 255.0d0 / 
     &              (maxTopChg-minTopChg)

               write(11,'(a)')    '***********************'
               write(11,'(a,3i3)') '3-Q separation scheme ', iiy, Tee, offset
               write(11,'(a)')    '***********************'
               write(11,'(f15.8)') actionA(:,:,:,iiy,Tee,offset)

               write(12,'(a)')    '***********************'
               write(12,'(a,3i3)') '3-Q separation scheme', iiy, Tee,offset
               write(12,'(a)')    '***********************'
               write(12,'(f15.8)') topChgA(:,:,:,iiy,Tee,offset)

            enddo               ! end iiy loop

         enddo                  ! end offset loop
      enddo                     ! end Tee loop

      close (11)
      close (12)
c
      close (9)
c
      end program AverageLT
