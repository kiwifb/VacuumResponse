c
c  A program to average the L loops results from VacuumRespLTplan
c  or VacuumRespLplan.
c  Adapted from AverageYall
c
c
      program AverageLall
c
      USE L_baryonParam
c
      IMPLICIT NONE
c
      integer,parameter                                       :: offmax=4
      integer                                                 :: start, finish, icon, Correlate
      integer                                                 :: iiy,offset,Tee
      integer                                                 :: ix,iy,iz
      integer                                                 :: icx,icy,icz
      integer                                                 :: ixq,iyq
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
c
      double precision,dimension(nL(that),offmax)             :: avgAction, avgTopChg
!HPF$ DISTRIBUTE avgAction(BLOCK,BLOCK)
!HPF$ DISTRIBUTE avgTopChg(BLOCK,BLOCK)
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
         write(*,'(/,3a,/)') 'Reading ',dialog,'correlations.'
c
         open (11,file=trim(thisConfig)//'.L1-xy'//trim(fstr1)//'.correl.unf',status='old',form='unformatted')
         open (12,file=trim(thisConfig)//'.L1-xy'//trim(fstr2)//'.correl.unf',status='old',form='unformatted')

         open (13,file=trim(thisConfig)//'.L2-xy'//trim(fstr1)//'.correl.unf',status='old',form='unformatted')
         open (14,file=trim(thisConfig)//'.L2-xy'//trim(fstr2)//'.correl.unf',status='old',form='unformatted')

         open (15,file=trim(thisConfig)//'.L3-xy'//trim(fstr1)//'.correl.unf',status='old',form='unformatted')
         open (16,file=trim(thisConfig)//'.L3-xy'//trim(fstr2)//'.correl.unf',status='old',form='unformatted')

         open (17,file=trim(thisConfig)//'.L4-xy'//trim(fstr1)//'.correl.unf',status='old',form='unformatted')
         open (18,file=trim(thisConfig)//'.L4-xy'//trim(fstr2)//'.correl.unf',status='old',form='unformatted')

         open (21,file=trim(thisConfig)//'.L1-xz'//trim(fstr1)//'.correl.unf',status='old',form='unformatted')
         open (22,file=trim(thisConfig)//'.L1-xz'//trim(fstr2)//'.correl.unf',status='old',form='unformatted')

         open (23,file=trim(thisConfig)//'.L2-xz'//trim(fstr1)//'.correl.unf',status='old',form='unformatted')
         open (24,file=trim(thisConfig)//'.L2-xz'//trim(fstr2)//'.correl.unf',status='old',form='unformatted')

         open (25,file=trim(thisConfig)//'.L3-xz'//trim(fstr1)//'.correl.unf',status='old',form='unformatted')
         open (26,file=trim(thisConfig)//'.L3-xz'//trim(fstr2)//'.correl.unf',status='old',form='unformatted')

         open (27,file=trim(thisConfig)//'.L4-xz'//trim(fstr1)//'.correl.unf',status='old',form='unformatted')
         open (28,file=trim(thisConfig)//'.L4-xz'//trim(fstr2)//'.correl.unf',status='old',form='unformatted')
c
c  ----------
c  accumulate
c  ----------
c
c  .L1-xy
c
         do ixq = 0, nL(xhat)
            do iyq = 0, nL(yhat)
               read(11) actionC(ixq,iyq,:,:,:,:,:)
               read(11) Wavg(ixq,iyq,:), avgAction(:,:)
               read(12) topChgC(ixq,iyq,:,:,:,:,:)
               read(12) Wavg(ixq,iyq,:), avgTopChg(:,:)
            end do
         end do
c
         ActionA    = ActionA    + ActionC
         topChgA    = topChgA    + topChgC
         avgActionA = avgActionA + avgAction
         avgTopChgA = avgTopChgA + avgTopChg
         WavgA      = WavgA      + Wavg
c
c  .L2-xy
c
         do ixq = 0, nL(xhat)
            do iyq = 0, nL(yhat)
               read(13) actionC(ixq,iyq,:,:,:,:,:)
               read(13) Wavg(ixq,iyq,:), avgAction(:,:)
               read(14) topChgC(ixq,iyq,:,:,:,:,:)
               read(14) Wavg(ixq,iyq,:), avgTopChg(:,:)
            end do
         end do
c
c  no need to add contribution from the L loop case where iyq=0
c  similarly, no need to add avgActionA and avgTopChgA again, but WavgA should.
c
c  however, care must be taken in the average step
c
         do iy = 0, ny-1

            icy = modulo(ny - iy , ny) 

            ActionA(:,1:nL(yhat),:,iy,:,:,:) = ActionA(:,1:nL(yhat),:, iy,:,:,:) + 
     &                                         ActionC(:,1:nl(yhat),:,icy,:,:,:)

            topChgA(:,1:nL(yhat),:,iy,:,:,:) = topChgA(:,1:nL(yhat),:, iy,:,:,:) + 
     &                                         topChgC(:,1:nL(yhat),:,icy,:,:,:)

         enddo
c
         WavgA(:,1:nL(yhat),:) = WavgA(:,1:nL(yhat),:) + Wavg(:,1:nL(yhat),:)
c
c  .L3-xy
c
         do ixq = 0, nL(xhat)
            do iyq = 0, nL(yhat)
               read(15) actionC(ixq,iyq,:,:,:,:,:)
               read(15) Wavg(ixq,iyq,:), avgAction(:,:)
               read(16) topChgC(ixq,iyq,:,:,:,:,:)
               read(16) Wavg(ixq,iyq,:), avgTopChg(:,:)
            end do
         end do
c
c  Here the only point that we don't need to add is the one for ixq=iyq=0
c
         do ix = 0,
c
c  .L2-xy
c
         do ixq = 0, nL(xhat)
            do iyq = 0, nL(yhat)
               read(17) actionC(ixq,iyq,:,:,:,:,:)
               read(17) Wavg(ixq,iyq,:), avgAction(:,:)
               read(18) topChgC(ixq,iyq,:,:,:,:,:)
               read(18) Wavg(ixq,iyq,:), avgTopChg(:,:)
            end do
         end do
c
c  no need to add contribution from the L loop case where iyq=0
c  similarly, no need to add avgActionA and avgTopChgA again, but WavgA should.
c
c  however, care must be taken in the average step
c
         do iy = 0, ny-1

            icy = modulo(ny - iy , ny) 

            ActionA(:,1:nL(yhat),:,iy,:,:,:) = ActionA(:,1:nL(yhat),:, iy,:,:,:) + 
     &                                         ActionC(:,1:nl(yhat),:,icy,:,:,:)

            topChgA(:,1:nL(yhat),:,iy,:,:,:) = topChgA(:,1:nL(yhat),:, iy,:,:,:) + 
     &                                         topChgC(:,1:nL(yhat),:,icy,:,:,:)

         enddo
c
         WavgA(:,1:nL(yhat),:) = WavgA(:,1:nL(yhat),:) + Wavg(:,1:nL(yhat),:)
c
c  .Tp-xz 
c
         do ixq = 0, nL(xhat)
            do iyq = 0, nL(yhat)
               read(21) actionC(ixq,iyq,:,:,:,:,:)
               read(21) Wavg(ixq,iyq,:), avgAction(:,:)
               read(22) topChgC(ixq,iyq,:,:,:,:,:)
               read(22) Wavg(ixq,iyq,:), avgTopChg(:,:)
            end do
         end do
c
         do iy = 0, nz-1
            do iz = 0, ny-1

               icy = modulo(nz - iz , nz)
               icz = modulo(     iy , ny)

               ActionA(:,1:nL(yhat),iz,iy,:,:,:) = ActionA(:,1:nL(yhat), iz, iy,:,:,:) + 
     &                                             ActionC(:,1:nL(yhat),icz,icy,:,:,:)

               topChgA(:,1:nL(yhat),iz,iy,:,:,:) = topChgA(:,1:nL(yhat), iz, iy,:,:,:) + 
     &                                             topChgC(:,1:nL(yhat),icz,icy,:,:,:)

            enddo
         enddo
c
         WavgA(:,1:nL(yhat),:) = WavgA(:,1:nL(yhat),:) + Wavg(:,1:nL(yhat),:)
c
c  .Tm-xz
c
         do ixq = 0, nL(xhat)
            do iyq = 0, nL(yhat)
               read(23) actionC(ixq,iyq,:,:,:,:,:)
               read(23) Wavg(ixq,iyq,:), avgAction(:,:)
               read(24) topChgC(ixq,iyq,:,:,:,:,:)
               read(24) Wavg(ixq,iyq,:), avgTopChg(:,:)
            end do
         end do
c
         do iy = 0, nz-1
            do iz = 0, ny-1
               do ix = 0, nt-1

                  icx = modulo(nt - ix , nt) 
                  icy = modulo(nz - iz , nz)
                  icz = modulo(     iy , ny)

                  ActionA(1:nL(xhat),1:nL(yhat),iz,iy,ix,:,:) = ActionA(1:nL(xhat),1:nL(yhat), iz, iy, ix,:,:) + 
     &                                                          ActionC(1:nL(xhat),1:nL(yhat),icz,icy,icx,:,:)

                  topChgA(1:nL(xhat),1:nL(yhat),iz,iy,ix,:,:) = topChgA(1:nL(xhat),1:nL(yhat), iz, iy, ix,:,:) + 
     &                                                          topChgC(1:nL(xhat),1:nL(yhat),icz,icy,icx,:,:)

               enddo
            enddo
         enddo
c
         WavgA(1:nL(xhat),1:nL(yhat),:) = WavgA(1:nL(xhat),1:nL(yhat),:) + Wavg(1:nL(xhat),1:nL(yhat),:)
c
         close(11)
         close(12)
         close(13)
         close(14)
         close(21)
         close(22)
         close(23)
         close(24)
c
      end do                    ! end icon loop
c
c  -------
c  average
c  -------
c
      actionA    = actionA    / (finish - start + 1)
      topChgA    = topChgA    / (finish - start + 1)
      WavgA      = WavgA      / (finish - start + 1)
      avgActionA = avgActionA / (finish - start + 1)
      avgTopChgA = avgTopChgA / (finish - start + 1)
c
c  We must also divide actionA, topChgA and Wavg by 4 for T loops indexed
c  by ixq=1,nL(xdir) and iyq=1,nL(ydir) 
c  since have added .Tp-xy, .Tm-xy, .Tp-xz and .Tm-xz
c
      actionA(1:nL(xhat),1:nL(yhat),:,:,:,:,:) = actionA(1:nL(xhat),1:nL(yhat),:,:,:,:,:) / 4
      topChgA(1:nL(xhat),1:nL(yhat),:,:,:,:,:) = topChgA(1:nL(xhat),1:nL(yhat),:,:,:,:,:) / 4
      WavgA(1:nL(xhat),1:nL(yhat),:) = WavgA(1:nL(xhat),1:nL(yhat),:) /4
c
c  We must also divide actionA, topChgA and Wavg by 2 for T loops indexed
c  by (ixq=0,iyq=1,nL(yhat)) and (ixq=1,nL(xhat),iyq=0) 
c  since have added .Tp-xy, .Tm-xy, .Tp-xz and .Tm-xz and we want to avoid double counting
c  In the first case (ixq=0) contributions from Tm and Tp are identical,
c  in the second (iyq=0) contributions T*-xy and T*-xz are identical.
c
      actionA(0,1:nL(yhat),:,:,:,:,:) = actionA(0,1:nL(yhat),:,:,:,:,:) / 2
      topChgA(0,1:nL(yhat),:,:,:,:,:) = topChgA(0,1:nL(yhat),:,:,:,:,:) / 2
      WavgA(0,1:nL(yhat),:) = WavgA(0,1:nL(yhat),:) /2
c
      actionA(1:nL(xhat),0,:,:,:,:,:) = actionA(1:nL(xhat),0,:,:,:,:,:) / 2
      topChgA(1:nL(xhat),0,:,:,:,:,:) = topChgA(1:nL(xhat),0,:,:,:,:,:) / 2
      WavgA(1:nL(xhat),0,:) = WavgA(1:nL(xhat),0,:) /2
c
c  And of course ixq=iyq=0 is the same in all four case and should only be counted once.
c  So there is no need to divide it by anything.
c
c  --------------------------------
c  move origin to center of lattice
c  --------------------------------
c
      actionA(:,:,:,:,:,:,:) = cshift(
     &     cshift(
     &     cshift(actionA(:,:,:,:,:,:,:), shift=-nt/2, dim=5),
     &     shift=-nz/2, dim=4),
     &     shift=-ny/2, dim=3)

      topChgA(:,:,:,:,:,:,:) = cshift(
     &     cshift(
     &     cshift(topChgA(:,:,:,:,:,:,:), shift=-nt/2, dim=5),
     &     shift=-nz/2, dim=4),
     &     shift=-ny/2, dim=3)
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
             
c     maxActionPref = maxval(actionA(:,:,:,3,Tee,offset))
c     minActionPref = minval(actionA(:,:,:,3,Tee,offset))
c     maxTopChgPref = maxval(topChgA(:,:,:,3,Tee,offset))
c     minTopChgPref = minval(topChgA(:,:,:,3,Tee,offset))

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
      end program AverageTall
