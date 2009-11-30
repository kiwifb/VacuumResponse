c
c  A program to average the results from VacuumRespLT
c  Specialized in getting the results from a T-shape
c  By F. Bissey Jan. 2004
c  Can also be used in conjunction with AverageTall to do
c  averaging in installement. In this case one has to make
c  sure that all partial average have the same weight, ie
c  they included the same number of configurations. 
c  centering the data is made optional as in the last case
c  the data may already be centered.
c
      program AverageT
c
      USE L_baryonParam
c
      IMPLICIT NONE
c
      integer,parameter                                       :: offmax=4
      integer                                                 :: start, finish, icon, Correlate, flag
      integer                                                 :: ixq,iyq,offset,Tee, center
      character(len=80)                                       :: configName, suffixName, reportName
      character(len=132)                                      :: thisConfig
      character(len=12)                                       :: fstr1,fstr2,wstr1,wstr2
      character(len=30)                                       :: dialog
c
      double precision                                        :: maxAction,  maxTopChg
      double precision                                        :: minAction,  minTopChg
      double precision                                        :: maxActionPref,  maxTopChgPref
      double precision                                        :: minActionPref,  minTopChgPref

      double precision,dimension(nL(that),offmax)             :: avgAction, avgTopChg, dummy
!HPF$ DISTRIBUTE avgAction(BLOCK,BLOCK)
!HPF$ DISTRIBUTE avgTopChg(BLOCK,BLOCK)
!HPF$ DISTRIBUTE dummy(BLOCK,BLOCK)
      double precision,dimension(nL(that),offmax)             :: avgActionA, avgTopChgA
!HPF$ DISTRIBUTE avgActionA(BLOCK,BLOCK)
!HPF$ DISTRIBUTE avgTopChgA(BLOCK,BLOCK)
c
      double precision,dimension(0:nL(xhat),0:nL(yhat),0:ny-1,0:nz-1,0:nt-1,nL(that),offmax)   :: actionC
!HPF$ DISTRIBUTE actionC(BLOCK,BLOCK,*,*,*,*,*)
      double precision,dimension(0:nL(xhat),0:nL(yhat),0:ny-1,0:nz-1,0:nt-1,nL(that),offmax)   :: topChgC
!HPF$ DISTRIBUTE topChgC(BLOCK,BLOCK,*,*,*,*,*)
      double precision,dimension(0:nL(xhat),0:nL(yhat),0:ny-1,0:nz-1,0:nt-1,nL(that),offmax)   :: actionA
!HPF$ DISTRIBUTE actionA(BLOCK,BLOCK,*,*,*,*,*)
      double precision,dimension(0:nL(xhat),0:nL(yhat),0:ny-1,0:nz-1,0:nt-1,nL(that),offmax)   :: topChgA
!HPF$ DISTRIBUTE topChgA(BLOCK,BLOCK,*,*,*,*,*)
      double precision,dimension(0:nL(xhat),0:nL(yhat),nL(that))  :: Wavg
!HPF$ DISTRIBUTE Wavg(BLOCK,BLOCK,*)
      double precision,dimension(0:nL(xhat),0:nL(yhat),nL(that))  :: WavgA
!HPF$ DISTRIBUTE WavgA(BLOCK,BLOCK,*)
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
c  The suffix is of the form ".s04" for example, it must include Txx if have it. 
c
      write(*,'(a)') 'What is the file name suffix?'
      read (*,'(a80)') suffixName
c
      write(*,'(a)') 'Please enter a report file name.'
      read (*,'(a80)') reportName
c
      write(*,'(a)') 'Do you want to center the data (0:yes, anything else no)?'
      read (*,*) center
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
c
c  Initialize accumulators
c
      actionA      = 0.0d0
      topChgA      = 0.0d0
      WavgA        = 0.0d0
      avgActionA   = 0.0d0
      avgTopChgA   = 0.0d0

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
         open (11,file=trim(thisConfig)//trim(fstr1)//'.correl.unf',status='old',form='unformatted')
         open (12,file=trim(thisConfig)//trim(fstr2)//'.correl.unf',status='old',form='unformatted')
c
c  ----------
c  accumulate
c  ----------
c
         do ixq = 0, nL(xhat)
            do iyq = 0, nL(yhat)
               do Tee = 2,nL(that)
                  do offset = 1,offmax

                     if ( Tee - 2*offset .lt. 0 ) then
                        cycle
                     endif

                     read(11) actionC(ixq,iyq,:,:,:,Tee,offset)
                     read(11) Wavg(ixq,iyq,Tee), avgAction(Tee,offset)
                     read(12) topChgC(ixq,iyq,:,:,:,Tee,offset)
                     read(12) Wavg(ixq,iyq,Tee), avgTopChg(Tee,offset)

                  end do
               end do
            end do
         end do
c
         ActionA    = ActionA    + ActionC
         topChgA    = topChgA    + topChgC
         avgActionA = avgActionA + avgAction
         avgTopChgA = avgTopChgA + avgTopChg
         WavgA      = WavgA      + Wavg

      end do
c
c  Average
c
      actionA      = actionA      / (finish - start + 1)
      topChgA      = topChgA      / (finish - start + 1)
      WavgA        = WavgA        / (finish - start + 1)
      avgActionA   = avgActionA   / (finish - start + 1)
      avgTopChgA   = avgTopChgA   / (finish - start + 1)

      if( center == 0 ) then
c
c  --------------------------------
c  move origin to center of lattice
c  --------------------------------
c
         actionA(:,:,:,:,:,:,:) = cshift(cshift(cshift(
     &        actionA(:,:,:,:,:,:,:), 
     &        shift=-nt/2, dim=5),
     &        shift=-nz/2, dim=4),
     &        shift=-ny/2, dim=3)

         topChgA(:,:,:,:,:,:,:) = cshift(cshift(cshift(
     &        topChgA(:,:,:,:,:,:,:), 
     &        shift=-nt/2, dim=5),
     &        shift=-nz/2, dim=4),
     &        shift=-ny/2, dim=3)
      endif
c
c  --------------
c  record results
c  --------------
c
      open ( 9,file=trim(reportName)//'.T.report',status='unknown')
c
      open (11,file=trim(reportName)//'.T'//trim(fstr1)//'.correl.unf',status='unknown',form='unformatted')
      open (12,file=trim(reportName)//'.T'//trim(fstr2)//'.correl.unf',status='unknown',form='unformatted')
c

         do ixq = 0, nL(xhat)
            do iyq = 0, nL(yhat)
               do Tee = 2,nL(that)
                  do offset = 1,offmax
c
                     if ( Tee - 2*offset .lt. 0 ) then
                        cycle
                     endif
c
                     write(11) actionA(ixq,iyq,:,:,:,Tee,offset)
                     write(11) WavgA(ixq,iyq,Tee), avgActionA(Tee,offset)
                     write(12) topChgA(ixq,iyq,:,:,:,Tee,offset)
                     write(12) WavgA(ixq,iyq,Tee), avgTopChgA(Tee,offset)
c
                  enddo                  ! end offset loop
               enddo                     ! end Tee loop
            end do
         end do
c
      close(11)
      close(12)
c
c  In AVS format unscaled
c
      open (11,file=trim(reportName)//'.T'//trim(fstr1)//'.correl',status='unknown')
      open (12,file=trim(reportName)//'.T'//trim(fstr2)//'.correl',status='unknown')
c
      do Tee = 2,nL(that)
         do offset = 1,offmax
            
            if ( Tee - 2*offset .lt. 0 ) then
               cycle
            endif
            
            do ixq = 0, nL(xhat)
               do iyq = 0, nL(yhat)

                  write(11,'(a)')    '*********************'
                  write(11,'(a,4i3)') 'Wilson Loop ', ixq, iyq, Tee,offset
                  write(11,'(a)')    '*********************'
                  write(11,'(f15.8)') actionA(ixq,iyq,:,:,:,Tee,offset) / (WavgA(ixq,iyq,Tee) * avgActionA(Tee,offset))

                  write(12,'(a)')    '*********************'
                  write(12,'(a,4i3)') 'Wilson Loop ', ixq, iyq, Tee,offset
                  write(12,'(a)')    '*********************'
                  write(12,'(f15.8)') topChgA(ixq,iyq,:,:,:,Tee,offset) / (WavgA(ixq,iyq,Tee) * avgTopChgA(Tee,offset))
            
               enddo            ! end iyq loop
            enddo               ! end ixq loop
            
         enddo                  ! end offset loop
      enddo                     ! end Tee loop

      close (11)
      close (12)
c
c  Normalized to 0 ... 255
c
      open (11,file=trim(reportName)//'.T'//trim(fstr1)//'.correl.255',status='unknown')
      open (12,file=trim(reportName)//'.T'//trim(fstr2)//'.correl.255',status='unknown')
c
c  Trap preferred values
c
      do Tee = 2,nL(that)
         do offset = 1,offmax

            if ( Tee - 2*offset .lt. 0 ) then
               cycle
            endif
             
            maxActionPref = 1.0082
            minActionPref = 0.9118
            maxTopChgPref = 1.0174
            minTopChgPref = 0.8478

            do ixq = 0, nL(xhat)
               do iyq = 0, nL(yhat) 

                  actionA(ixq,iyq,:,:,:,Tee,offset) = 
     &                 actionA(ixq,iyq,:,:,:,Tee,offset) / (WavgA(ixq,iyq,Tee) * avgActionA(Tee,offset))
                   
                  topChgA(ixq,iyq,:,:,:,Tee,offset) = 
     &                 topChgA(ixq,iyq,:,:,:,Tee,offset) / (WavgA(ixq,iyq,Tee) * avgTopChgA(Tee,offset))
c
c  report max min values of correlations
c
                  write(9,'(/,a,4i3,/,4(3a,f8.4,/))') 
     &                 'At separation, ', ixq,iyq,Tee,offset,
     &                 'Max',wstr1,'Correlation = ', maxval(actionA(ixq,iyq,:,:,:,Tee,offset)),
     &                 'Min',wstr1,'Correlation = ', minval(actionA(ixq,iyq,:,:,:,Tee,offset)),
     &                 'Max',wstr2,'Correlation = ', maxval(topChgA(ixq,iyq,:,:,:,Tee,offset)),
     &                 'Min',wstr2,'Correlation = ', minval(topChgA(ixq,iyq,:,:,:,Tee,offset))

                  maxAction = max(maxval(actionA(ixq,iyq,:,:,:,Tee,offset)),maxActionPref)
                  minAction = min(minval(actionA(ixq,iyq,:,:,:,Tee,offset)),minActionPref)
                  actionA(ixq,iyq,:,:,:,Tee,offset) = (actionA(ixq,iyq,:,:,:,Tee,offset)-minAction) * 255.0d0 / 
     &                 (maxAction - minAction)

                  maxTopChg = max(maxval(topChgA(ixq,iyq,:,:,:,Tee,offset)),maxTopChgPref)
                  minTopChg = min(minval(topChgA(ixq,iyq,:,:,:,Tee,offset)),minTopChgPref)
                  topChgA(ixq,iyq,:,:,:,Tee,offset) = (topChgA(ixq,iyq,:,:,:,Tee,offset)-minTopChg) * 255.0d0 / 
     &                 (maxTopChg-minTopChg)

                  write(11,'(a)')    '*********************'
                  write(11,'(a,4i3)') 'Q Q-bar separation ', ixq,iyq, Tee, offset
                  write(11,'(a)')    '*********************'
                  write(11,'(f15.8)') actionA(ixq,iyq,:,:,:,Tee,offset)

                  write(12,'(a)')    '*********************'
                  write(12,'(a,4i3)') 'Q Q-bar separation ', ixq,iyq, Tee,offset
                  write(12,'(a)')    '*********************'
                  write(12,'(f15.8)') topChgA(ixq,iyq,:,:,:,Tee,offset)

               enddo            ! end iyq loop
            enddo               ! end ixq loop

         enddo                  ! end offset loop
      enddo                     ! end Tee loop

      close (11)
      close (12)
c
      close (9)


c
      end program AverageT
