c
c  A program to average the results from VacuumRespLT
c  Specialized in getting the results from the T-shape
c  By F. Bissey Jan. 2004
c
      program AverageTrd
c
      USE L_baryonParam
c
      IMPLICIT NONE
c
      integer,parameter                                       :: offmax=4
      integer                                                 :: start, finish, icon, Correlate, flag
      integer                                                 :: ixq,iyq,offset,Tee
      character(len=3)                                        :: cfg
      character(len=80)                                       :: configName, suffixName, reportName
      character(len=132)                                      :: thisConfig
      character(len=12)                                       :: fstr1,fstr2,wstr1,wstr2
      character(len=30)                                       :: dialog,loop
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

      double precision,dimension(0:ny-1,0:nz-1,0:nt-1,nL(that),offmax)   :: actionC
!HPF$ DISTRIBUTE actionC(*,BLOCK,BLOCK,*,*)
      double precision,dimension(0:ny-1,0:nz-1,0:nt-1,nL(that),offmax)   :: topChgC
!HPF$ DISTRIBUTE topChgC(*,BLOCK,BLOCK,*,*)
      double precision,dimension(0:ny-1,0:nz-1,0:nt-1,-nL(xhat):nL(xhat),-nL(yhat):nL(yhat),nL(that),offmax)   :: actionA
!HPF$ DISTRIBUTE actionA(*,BLOCK,BLOCK,*,*,*,*)
      double precision,dimension(0:ny-1,0:nz-1,0:nt-1,-nL(xhat):nL(xhat),-nL(yhat):nL(yhat),nL(that),offmax)   :: topChgA
!HPF$ DISTRIBUTE topChgA(*,BLOCK,BLOCK,*,*,*,*)
      double precision,dimension(nL(that))                     :: Wavg
!HPF$ DISTRIBUTE Wavg(*)
      double precision,dimension(-nL(xhat):nL(xhat),-nL(yhat):nL(yhat),nL(that))   :: WavgA
!HPF$ DISTRIBUTE WavgA(*,*,*)

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
c T or L shape loop? The "p" and "m" will be added appropriately    
c
      write(*,'(a)') 'What kind of loop? [T or L, in caps please]'
      read(*,'(a30)') loop
 
      write(*,'(a)') 'Please enter a report file name.'
      read (*,'(a80)') reportName
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
c=============================================================================
c     Compute p and m
c=============================================================================
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
         write(thisConfig,fmt='(a,i3.3,a,a,a)') trim(configName),icon,trim(suffixName),'.',trim(loop)
         
         write(*,*) trim(thisConfig)
c
c  Read in results for one configuration
c
         write(*,'(/,3a,/)') 'Reading ',dialog,'correlations.'
         open (11,file=trim(thisConfig)//'p'//trim(fstr1)//'.correl.unf',status='old',form='unformatted')
         open (12,file=trim(thisConfig)//'p'//trim(fstr2)//'.correl.unf',status='old',form='unformatted')
         open (13,file=trim(thisConfig)//'m'//trim(fstr1)//'.correl.unf',status='old',form='unformatted')
         open (14,file=trim(thisConfig)//'m'//trim(fstr2)//'.correl.unf',status='old',form='unformatted')
c
c   avgAction and avgTopChg are recorded for each value of ixq and iyq
c   however they are independent from it. Rather than dividing the final
c   values of avgActionA and avgTopChgA by the number of time they are counted
c   we prefer to accumulate them only the first time they are encountered 
c   in the record.  
c==============================================================================
c
c   For ixq=0 the Tpp and Tmm shapes are the same. So only one of them should
c   added to the average.
c
         flag = 0

         do ixq = 0, nL(xhat)
            do iyq= 0, nL(yhat)
               
               read(11) actionC(:,:,:,:,:)
               read(11) Wavg(:), avgAction(:,:)
               
               read(12) topChgC(:,:,:,:,:)
               read(12) Wavg(:), avgTopChg(:,:)
c
c Accumulate
c
               actionA(:,:,:,ixq,iyq,:,:) = actionA(:,:,:,ixq,iyq,:,:) + actionC(:,:,:,:,:)
               topChgA(:,:,:,ixq,iyq,:,:) = topChgA(:,:,:,ixq,iyq,:,:) + topChgC(:,:,:,:,:)

               WavgA(ixq,iyq,:) = WavgA(ixq,iyq,:) + Wavg(:)

               if ( flag == 0) then 

                  avgActionA = avgActionA + avgAction
                  avgTopChgA = avgTopChgA + avgTopChg
                  
               endif

               read(13) actionC(:,:,:,:,:)
               read(13) Wavg(:), avgAction(:,:)
               
               read(14) topChgC(:,:,:,:,:)
               read(14) Wavg(:), avgTopChg(:,:)
c
c Accumulate
c
               if ( ixq .ne. 0 ) then

                  actionA(:,:,:,-ixq,iyq,:,:) = actionA(:,:,:,-ixq,iyq,:,:) + actionC(:,:,:,:,:)
                  topChgA(:,:,:,-ixq,iyq,:,:) = topChgA(:,:,:,-ixq,iyq,:,:) + topChgC(:,:,:,:,:)

                  WavgA(-ixq,iyq,:) = WavgA(-ixq,iyq,:) + Wavg(:)

               end if
               
               flag = 1         ! set flag to 1 for all subsequent passage in the heart of the loop

            enddo               ! end iyq loop
         enddo                  ! end ixq loop
         
         close(11)
         close(12)
         close(13)
         close(14)
           
      end do                    ! end icon loop
c
c  Average
c
      actionA      = actionA      / (finish - start + 1)
      topChgA      = topChgA      / (finish - start + 1)
      WavgA        = WavgA        / (finish - start + 1)
      avgActionA   = avgActionA   / (finish - start + 1)
      avgTopChgA   = avgTopChgA   / (finish - start + 1)

      actionA(:,:,:,:,:,:,:) = cshift(
     &     cshift(
     &     cshift(actionA(:,:,:,:,:,:,:), shift=-nt/2, dim=3),
     &     shift=-nz/2, dim=2),
     &     shift=-ny/2, dim=1)

      topChgA(:,:,:,:,:,:,:) = cshift(
     &     cshift(
     &     cshift(topChgA(:,:,:,:,:,:,:), shift=-nt/2, dim=3),
     &     shift=-nz/2, dim=2),
     &     shift=-ny/2, dim=1)

      open (9,file=trim(reportName)//'.'//trim(loop)//'.report',status='unknown')

      open (11,file=trim(reportName)//'.'//trim(loop)//trim(fstr1)//'.correl.unf',status='unknown',form='unformatted')
      open (12,file=trim(reportName)//'.'//trim(loop)//trim(fstr2)//'.correl.unf',status='unknown',form='unformatted')

      do Tee = 2,nL(that)
         do offset = 1,offmax
 
            if ( Tee - 2*offset .lt. 0 ) then
               cycle
            endif

            write(11) actionA(:,:,:,:,0:nL(yhat),Tee,offset)
            write(11) WavgA(:,0:nL(yhat),Tee), avgActionA(Tee,offset)
            write(12) topChgA(:,:,:,:,0:nL(yhat),Tee,offset)
            write(12) WavgA(:,0:nL(yhat),Tee), avgTopChgA(Tee,offset)

         enddo                  ! end offset loop
      enddo                     ! end Tee loop

      close(11)
      close(12)
c
c  In AVS format unscaled
c
      open (11,file=trim(reportName)//'.'//trim(loop)//trim(fstr1)//'.correl',status='unknown')
      open (12,file=trim(reportName)//'.'//trim(loop)//trim(fstr2)//'.correl',status='unknown')
c
      do Tee = 2,nL(that)
         do offset = 1,offmax
            
            if ( Tee - 2*offset .lt. 0 ) then
               cycle
            endif
            
            do iyq = 0, nL(yhat)
               do ixq = -6, 6

                  write(11,'(a)')    '*********************'
                  write(11,'(a,4i3)') 'Wilson Loop ', ixq, iyq, Tee,offset
                  write(11,'(a)')    '*********************'
                  write(11,'(f15.8)') actionA(:,:,:,ixq,iyq,Tee,offset) / (WavgA(ixq,iyq,Tee) * avgActionA(Tee,offset))

                  write(12,'(a)')    '*********************'
                  write(12,'(a,4i3)') 'Wilson Loop ', ixq, iyq, Tee,offset
                  write(12,'(a)')    '*********************'
                  write(12,'(f15.8)') topChgA(:,:,:,ixq,iyq,Tee,offset) / (WavgA(ixq,iyq,Tee) * avgTopChgA(Tee,offset))
            
               enddo            ! end ixq loop
            enddo               ! end iyq loop
            
         enddo                  ! end offset loop
      enddo                     ! end Tee loop

      close (11)
      close (12)
c
c  Normalized to 0 ... 255
c
      open (11,file=trim(reportName)//'.'//trim(loop)//trim(fstr1)//'.correl.255',status='unknown')
      open (12,file=trim(reportName)//'.'//trim(loop)//trim(fstr2)//'.correl.255',status='unknown')
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

            do iyq = 0, nL(yhat)
               do ixq = -6, 6 

                  actionA(:,:,:,ixq,iyq,Tee,offset) = 
     &                 actionA(:,:,:,ixq,iyq,Tee,offset) / (WavgA(ixq,iyq,Tee) * avgActionA(Tee,offset))
                   
                  topChgA(:,:,:,ixq,iyq,Tee,offset) = 
     &                 topChgA(:,:,:,ixq,iyq,Tee,offset) / (WavgA(ixq,iyq,Tee) * avgTopChgA(Tee,offset))
c
c  report max min values of correlations
c
                  write(9,'(/,a,4i3,/,4(3a,f8.4,/))') 
     &                 'At separation, ', ixq,iyq,Tee,offset,
     &                 'Max',wstr1,'Correlation = ', maxval(actionA(:,:,:,ixq,iyq,Tee,offset)),
     &                 'Min',wstr1,'Correlation = ', minval(actionA(:,:,:,ixq,iyq,Tee,offset)),
     &                 'Max',wstr2,'Correlation = ', maxval(topChgA(:,:,:,ixq,iyq,Tee,offset)),
     &                 'Min',wstr2,'Correlation = ', minval(topChgA(:,:,:,ixq,iyq,Tee,offset))

                  maxAction = max(maxval(actionA(:,:,:,ixq,iyq,Tee,offset)),maxActionPref)
                  minAction = min(minval(actionA(:,:,:,ixq,iyq,Tee,offset)),minActionPref)
                  actionA(:,:,:,ixq,iyq,Tee,offset) = (actionA(:,:,:,ixq,iyq,Tee,offset)-minAction) * 255.0d0 / 
     &                 (maxAction - minAction)

                  maxTopChg = max(maxval(topChgA(:,:,:,ixq,iyq,Tee,offset)),maxTopChgPref)
                  minTopChg = min(minval(topChgA(:,:,:,ixq,iyq,Tee,offset)),minTopChgPref)
                  topChgA(:,:,:,ixq,iyq,Tee,offset) = (topChgA(:,:,:,ixq,iyq,Tee,offset)-minTopChg) * 255.0d0 / 
     &                 (maxTopChg-minTopChg)

                  write(11,'(a)')    '*********************'
                  write(11,'(a,4i3)') 'Q Q-bar separation ', ixq,iyq, Tee, offset
                  write(11,'(a)')    '*********************'
                  write(11,'(f15.8)') actionA(:,:,:,ixq,iyq,Tee,offset)

                  write(12,'(a)')    '*********************'
                  write(12,'(a,4i3)') 'Q Q-bar separation ', ixq,iyq, Tee,offset
                  write(12,'(a)')    '*********************'
                  write(12,'(f15.8)') topChgA(:,:,:,ixq,iyq,Tee,offset)

               enddo            ! end ixq loop
            enddo               ! end iyq loop

         enddo                  ! end offset loop
      enddo                     ! end Tee loop

      close (11)
      close (12)
c
      close (9)
c
c=============================================================================
c-----------------------------------------------------------------------------
c=============================================================================

c
      end program AverageTrd
