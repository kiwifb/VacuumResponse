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
c
      program AverageYall
c
      USE L_baryonParam
      USE L_WRITEAVG
c
      IMPLICIT NONE
c
      include'Yfiles/yLoopSize.h'
      integer,parameter                                       :: offmax=4
      integer                                                 :: start, finish, icon, Correlate
      integer                                                 :: iiy,offset,Tee,icon2
      integer                                                 :: ix,iy,iz
      integer                                                 :: icx,icy,icz
      integer                                                 :: YorTY,iylp,av4,pot,zero
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
      double precision,dimension(0:ny-1,0:nz-1,0:nt-1,0:nYloop,nL(that),offmax)   :: actionC
!HPF$ DISTRIBUTE actionC(*,BLOCK,BLOCK,*,*,*)
      double precision,dimension(0:ny-1,0:nz-1,0:nt-1,0:nYloop,nL(that),offmax)   :: topChgC
!HPF$ DISTRIBUTE topChgC(*,BLOCK,BLOCK,*,*,*)
      double precision,dimension(0:ny-1,0:nz-1,0:nt-1,0:nYloop,nL(that),offmax)   :: actionA
!HPF$ DISTRIBUTE actionA(*,BLOCK,BLOCK,*,*,*)
      double precision,dimension(0:ny-1,0:nz-1,0:nt-1,0:nYloop,nL(that),offmax)   :: topChgA
!HPF$ DISTRIBUTE topChgA(*,BLOCK,BLOCK,*,*,*)
      double precision,dimension(0:nYloop,nL(that))           :: Wavg
!HPF$ DISTRIBUTE Wavg(BLOCK,BLOCK)
      double precision,dimension(0:nYloop,nL(that))           :: WavgA,WavgA2
!HPF$ DISTRIBUTE WavgA(BLOCK,BLOCK)
      integer,dimension(nYloop)                               :: shx
!HPF$ DISTRIBUTE shx(*)
      double precision,dimension(0:nYloop,nL(that))           :: WL
!HPF$ DISTRIBUTE WL(BLOCK,BLOCK)
      double precision,dimension(0:nYloop,nL(that),200)       :: wg,Vi
      double precision,dimension(0:nYloop,nL(that)-1)         :: V0,Vb,Vberr
      double precision,dimension(nYloop)                      :: rs, dqq
      integer,dimension(nYloop,nL(that))                      :: ncount
c
c  Begin execution
c
      rs(1) = 0.15301
      rs(2) = 0.26503
      rs(3) = 0.30603
      rs(4) = 0.41804
      rs(5) = 0.57106
      rs(6) = 0.72407
      rs(7) = 0.87708

      dqq(1) = 0.25568
      dqq(2) = 0.39383
      dqq(3) = 0.51136
      dqq(4) = 0.64907
      dqq(5) = 0.90455
      dqq(6) = 1.16012
      dqq(7) = 1.41573
      
      shx(1)=1
      shx(2)=1
      shx(3)=1
      shx(4)=2
      shx(5)=3
      shx(6)=4
      shx(7)=4
      
      if( nYloop .gt. 7) then
         
         rs(8) = 0.9891
         rs(9) = 1.14211
         
         dqq(8) = 1.55359
         dqq(9) = 1.80911
         
         shx(8)=4
         shx(9)=5
         
      end if      
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
      write(*,'(a)') 'Are we dealing with Y(1) or TY(2) results?'
      read (*,*) YorTY
c
      write(*,'(a)') 'Are we dealing Avg4 results (yes=1 no=0)?'
      read (*,*) av4
c
      write(*,'(a)') 'Do you want a potential study with that (yes=1 no=0)?'
      write(*,'(a)') '(Will not work if you do an averaging by part.)'
      read (*,*) pot
c
      write(*,'(a)') 'are averaging by part (yes=1 no=0)?'
      read (*,*) zero
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
      wg           = 0.0d0
      WavgA2       = 0.0d0
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
      if(YorTY==1) then
         write(suffixName,fmt='(a,a)') trim(suffixName),'.Y'
         if(zero == 0) then
            write(reportName,fmt='(a,a)') trim(reportName),'.Y'
         end if      
      else
         write(suffixName,fmt='(a,a)') trim(suffixName),'.TY'
         if(zero == 0) then
            write(reportName,fmt='(a,a)') trim(reportName),'.TY'
         end if
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
         if(av4 == 0) then
            
            open (11,file=trim(thisConfig)//'p-xy'//trim(fstr1)//'.correl.unf',status='old',form='unformatted')
            open (12,file=trim(thisConfig)//'p-xy'//trim(fstr2)//'.correl.unf',status='old',form='unformatted')
            open (13,file=trim(thisConfig)//'m-xy'//trim(fstr1)//'.correl.unf',status='old',form='unformatted')
            open (14,file=trim(thisConfig)//'m-xy'//trim(fstr2)//'.correl.unf',status='old',form='unformatted')

            open (21,file=trim(thisConfig)//'p-xz'//trim(fstr1)//'.correl.unf',status='old',form='unformatted')
            open (22,file=trim(thisConfig)//'p-xz'//trim(fstr2)//'.correl.unf',status='old',form='unformatted')
            open (23,file=trim(thisConfig)//'m-xz'//trim(fstr1)//'.correl.unf',status='old',form='unformatted')
            open (24,file=trim(thisConfig)//'m-xz'//trim(fstr2)//'.correl.unf',status='old',form='unformatted')
         
         else
         
            open (11,file=trim(thisConfig)//'avg4'//trim(fstr1)//'.correl.unf',status='old',form='unformatted')
            open (12,file=trim(thisConfig)//'avg4'//trim(fstr2)//'.correl.unf',status='old',form='unformatted')
         
         end if            
c
c  ----------
c  accumulate
c  ----------
c
c  .Yp-xy or avg4
c
         read(11) actionC(:,:,:,:,:,:)
         read(11) Wavg(:,:), avgAction(:,:)
         read(12) topChgC(:,:,:,:,:,:)
         read(12) Wavg(:,:), avgTopChg(:,:)

         if((YorTY == 2).and.(av4 == 0).and.(zero == 0)) then
            do iylp = 1,nYloop
         
               actionC(:,:,:,iylp,:,:) = cshift(actionC(:,:,:,iylp,:,:),dim=3,shift=-shx(iylp))

               topChgC(:,:,:,iylp,:,:) = cshift(topChgC(:,:,:,iylp,:,:),dim=3,shift=-shx(iylp))

            end do
         end if
c
         ActionA    = ActionA    + ActionC
         topChgA    = topChgA    + topChgC
         avgActionA = avgActionA + avgAction
         avgTopChgA = avgTopChgA + avgTopChg
         WavgA      = WavgA      + Wavg
         WL         = Wavg

         do iylp = 1, nYloop
            do Tee = 1, nL(that)

               if(Wavg(iylp,Tee)>0) then
                  WavgA2(iylp,Tee) = WavgA2(iylp,Tee) + Wavg(iylp,Tee)
                  ncount(iylp,Tee) = ncount(iylp,Tee) + 1
               end if

            end do
         end do

         if(av4 == 0) then
c
c  .Ym-xy
c
            read(13) actionC(:,:,:,:,:,:)
            read(13) Wavg(:,:), avgAction(:,:)
            read(14) topChgC(:,:,:,:,:,:)
            read(14) Wavg(:,:), avgTopChg(:,:)

            if(YorTY==2) then
               do iylp = 1,nYloop
         
                  actionC(:,:,:,iylp,:,:) = cshift(actionC(:,:,:,iylp,:,:),dim=3,shift=-shx(iylp))

                  topChgC(:,:,:,iylp,:,:) = cshift(topChgC(:,:,:,iylp,:,:),dim=3,shift=-shx(iylp)) 

               end do
            end if
c
c  no need to add contribution from first Y loop case (all three quarks at (0,0))
c  similarly, no need to add avgActionA and avgTopChgA again, but WavgA should.
c
c  however, care must be taken in the average step
c
            do ix = 0, nt-1

               icx = modulo(nt - ix , nt) 

               ActionA(:,:,ix,1:nYloop,:,:) = ActionA(:,:, ix,1:nYloop,:,:) + ActionC(:,:,icx,1:nYloop,:,:)

               topChgA(:,:,ix,1:nYloop,:,:) = topChgA(:,:, ix,1:nYloop,:,:) + topChgC(:,:,icx,1:nYloop,:,:)

            end do
c
            WavgA(1:nYloop,:) = WavgA(1:nYloop,:) + Wavg(1:nYloop,:)
            WL                = WL + Wavg
c
c  .Yp-xz 
c
            read(21) actionC(:,:,:,:,:,:)
            read(21) Wavg(:,:), avgAction(:,:)
            read(22) topChgC(:,:,:,:,:,:)
            read(22) Wavg(:,:), avgTopChg(:,:)

            if(YorTY==2) then
               do iylp = 1,nYloop
         
                  actionC(:,:,:,iylp,:,:) = cshift(actionC(:,:,:,iylp,:,:),dim=3,shift=-shx(iylp))
 
                  topChgC(:,:,:,iylp,:,:) = cshift(topChgC(:,:,:,iylp,:,:),dim=3,shift=-shx(iylp)) 

               end do
            end if
c
            do iy = 0, nz-1
               do iz = 0, ny-1

                  icy = modulo(nz - iz , nz)
                  icz = modulo(     iy , ny)

                  ActionA(iz,iy,:,1:nYloop,:,:) = ActionA(iz,iy,:,1:nYloop,:,:) + ActionC(icz,icy,:,1:nYloop,:,:)

                  topChgA(iz,iy,:,1:nYloop,:,:) = topChgA(iz,iy,:,1:nYloop,:,:) + topChgC(icz,icy,:,1:nYloop,:,:)

               end do
            end do
c
            WavgA(1:nYloop,:) = WavgA(1:nYloop,:) + Wavg(1:nYloop,:)
            WL                = WL + Wavg
c
c  .Ym-xz
c
            read(23) actionC(:,:,:,:,:,:)
            read(23) Wavg(:,:), avgAction(:,:)
            read(24) topChgC(:,:,:,:,:,:)
            read(24) Wavg(:,:), avgTopChg(:,:)

            if(YorTY==2) then
               do iylp = 1,nYloop
         
                  actionC(:,:,:,iylp,:,:) = cshift(actionC(:,:,:,iylp,:,:),dim=3,shift=-shx(iylp))

                  topChgC(:,:,:,iylp,:,:) = cshift(topChgC(:,:,:,iylp,:,:),dim=3,shift=-shx(iylp)) 
            
               end do
            end if
c
            do iy = 0, nz-1
               do iz = 0, ny-1
                  do ix = 0, nt-1

                     icx = modulo(nt - ix , nt) 
                     icy = modulo(nz - iz , nz)
                     icz = modulo(     iy , ny)

                     ActionA(iz,iy,ix,1:nYloop,:,:) = ActionA(iz,iy,ix,1:nYloop,:,:) + ActionC(icz,icy,icx,1:nYloop,:,:)

                     topChgA(iz,iy,ix,1:nYloop,:,:) = topChgA(iz,iy,ix,1:nYloop,:,:) + topChgC(icz,icy,icx,1:nYloop,:,:)

                  end do
               end do
            end do
c
            WavgA(1:nYloop,:) = WavgA(1:nYloop,:) + Wavg(1:nYloop,:)
            WL                = WL + Wavg
c
            WL = WL / 4
c
         end if
         
         close(11)
         close(12)
         
         if(av4 == 0) then
         
            close(13)
            close(14)
            close(21)
            close(22)
            close(23)
            close(24)
            
         end if
c
c----------------------------------------
c potential computations
c----------------------------------------
c
         if(pot == 1) then

            do icon2 = start, finish
            
               if(icon2 /= icon) then
                  wg(1:nYloop,:,icon2) = wg(1:nYloop,:,icon2) + WL(1:nYloop,:)
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
      actionA    = actionA    / (finish - start + 1 )
      topChgA    = topChgA    / (finish - start + 1 )
      WavgA      = WavgA      / (finish - start + 1 )
      avgActionA = avgActionA / (finish - start + 1 )
      avgTopChgA = avgTopChgA / (finish - start + 1 )


      WavgA2(1:nYloop,1:nL(that)) = WavgA2(1:nYloop,1:nL(that)) / ncount(1:nYloop,1:nL(that))

      if(pot == 1) then
         wg      = wg         / (finish - start )
      end if
      
      if(av4 == 0) then
c
c  must also divide actionA, topChgA and Wavg by 4 for Y loops 1 to nYloop 
c  since have added .Yp-xy, .Ym-xy, .Yp-xz and .Ym-xz
c
         actionA(:,:,:,1:nYloop,:,:) = actionA(:,:,:,1:nYloop,:,:) / 4
         topChgA(:,:,:,1:nYloop,:,:) = topChgA(:,:,:,1:nYloop,:,:) / 4
         WavgA(1:nYloop,:) = WavgA(1:nYloop,:) /4

      end if

      if(zero==0) then
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
c
c  ---------------------------------------
c  move to the loop center for comparison
c  if we have TY shapes.
c  ---------------------------------------
c
         if((YorTY == 2).and.(av4 == 0)) then
            do iylp = 1,nYloop
         
               actionA(:,:,:,iylp,:,:) = cshift(actionA(:,:,:,iylp,:,:),dim=3,shift=shx(iylp))

               topChgA(:,:,:,iylp,:,:) = cshift(topChgA(:,:,:,iylp,:,:),dim=3,shift=shx(iylp))

            end do
         endif

         open (9,file=trim(reportName)//'.wavg2',status='unknown')

         do iylp = 1, nYloop

            write(9,*) WavgA2(iylp,:)
            write(9,*) ncount(iylp,:)

         end do

         close(9)

c
c----------------------------------
c potential computations
c----------------------------------
c
         if(pot == 1) then

            do icon = start, finish
            
               do Tee = 1, nL(that) - 1

                  do iiy = 1, nYloop

                     Vi(iiy,Tee,icon) = log(wg(iiy,Tee,icon)/wg(iiy,Tee+1,icon))
                  
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

                  Vb(iiy,Tee) = sum(Vi(iiy,Tee,:)) /(finish - start +1)
                  Vberr(iiy,Tee) = sqrt(
     &                 ((sum((Vi(iiy,Tee,:)-Vb(iiy,Tee))**2))/
     &                 (finish - start +1))*(finish - start))

               end do
            end do   
      
            do Tee = 1,nL(that)
               do iiy = 1, nYloop

                  Vb(iiy,Tee) = log(WavgA2(iiy,Tee)/WavgA2(iiy,Tee+1))

               end do
            end do

         end if
c
c  --------------
c  record results
c  --------------
c
c----------------------------------
c record potential results if any
c----------------------------------
c
         if(pot == 1) then

            if(YorTY==1) then
               write(suffixName,fmt='(a)') '.wilson-Y.dat'
            else
               write(suffixName,fmt='(a)') '.wilson-T.dat'
            endif

            do iiy = 1, nYloop

               write(thisconfig,fmt='(a,a,i2.2,a)') trim(reportName),'-sh',iiy,trim(suffixName) 
        
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

            write(thisconfig,fmt='(a,a,a)') trim(reportName),'-pot',trim(suffixName) 
      
            open (11,file=trim(thisconfig),status='unknown',form='formatted')
      
            do iiy = 1, nYloop

               write(11,*) rs(iiy),dqq(iiy),V0(iiy,1),Vberr(iiy,1),Vb(iiy,1)

            end do
      
            close(11)
         
            write(thisconfig,fmt='(a,a,a)') trim(reportName),'-pot-T2',trim(suffixName) 
      
            open (11,file=trim(thisconfig),status='unknown',form='formatted')
      
            do iiy = 1, nYloop

               if(WavgA(iiy,3)<0.d0) then 
                  cycle
               end if

               write(11,*) rs(iiy),dqq(iiy),V0(iiy,2),Vberr(iiy,2),Vb(iiy,2)

            end do
      
            close(11)
         
            write(thisconfig,fmt='(a,a,a)') trim(reportName),'-pot-T3',trim(suffixName) 
      
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
c  should be the same as original AverageY?
c
         open (9,file=trim(reportName)//'.report',status='unknown')
c
         open (11,file=trim(reportName)//trim(fstr1)//'.correl.unf',status='unknown',form='unformatted')
         open (12,file=trim(reportName)//trim(fstr2)//'.correl.unf',status='unknown',form='unformatted')
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
         open (11,file=trim(reportName)//trim(fstr1)//'.correl',status='unknown')
         open (12,file=trim(reportName)//trim(fstr2)//'.correl',status='unknown')
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
                  write(11,'(a,3i3)') 'Y Loop ', iiy, Tee,offset
                  write(11,'(a)')    '*********************'
                  call fmt_writeavg(11,actionA(:,:,:,iiy,Tee,offset) / (WavgA(iiy,Tee) * avgActionA(Tee,offset)))
c
                  write(12,'(a)')    '*********************'
                  write(12,'(a,3i3)') 'Y Loop ', iiy, Tee,offset
                  write(12,'(a)')    '*********************'
                  call fmt_writeavg(12,topChgA(:,:,:,iiy,Tee,offset) / (WavgA(iiy,Tee) * avgTopChgA(Tee,offset)))
c           
               enddo               ! end iiy loop
c            
            enddo                  ! end offset loop
         enddo                     ! end Tee loop
c
         close (11)
         close (12)
c
c  Normalized to 0 ... 255
c
         open (11,file=trim(reportName)//trim(fstr1)//'.correl.255',status='unknown')
         open (12,file=trim(reportName)//trim(fstr2)//'.correl.255',status='unknown')
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
                  actionA(:,:,:,iiy,Tee,offset) = 
     &                 actionA(:,:,:,iiy,Tee,offset) / (WavgA(iiy,Tee) * avgActionA(Tee,offset))
c                   
                  topChgA(:,:,:,iiy,Tee,offset) = 
     &                 topChgA(:,:,:,iiy,Tee,offset) / (WavgA(iiy,Tee) * avgTopChgA(Tee,offset))
c
c  report max min values of correlations
c     for action and topological charge
                  write(9,'(/,a,3i3,/,4(3a,f8.4,/))') 
     &                 'At separation scheme, ', iiy,Tee,offset,
     &                 'Max',wstr1,'Correlation = ', maxval(actionA(:,:,:,iiy,Tee,offset)),
     &                 'Min',wstr1,'Correlation = ', minval(actionA(:,:,:,iiy,Tee,offset)),
     &                 'Max',wstr2,'Correlation = ', maxval(topChgA(:,:,:,iiy,Tee,offset)),
     &                 'Min',wstr2,'Correlation = ', minval(topChgA(:,:,:,iiy,Tee,offset))
c
                  maxAction = max(maxval(actionA(:,:,:,iiy,Tee,offset)),maxActionPref)
                  minAction = min(minval(actionA(:,:,:,iiy,Tee,offset)),minActionPref)
                  actionA(:,:,:,iiy,Tee,offset) = (actionA(:,:,:,iiy,Tee,offset)-minAction) * 255.0d0 / 
     &                 (maxAction - minAction)
c
                  maxTopChg = max(maxval(topChgA(:,:,:,iiy,Tee,offset)),maxTopChgPref)
                  minTopChg = min(minval(topChgA(:,:,:,iiy,Tee,offset)),minTopChgPref)
                  topChgA(:,:,:,iiy,Tee,offset) = (topChgA(:,:,:,iiy,Tee,offset)-minTopChg) * 255.0d0 / 
     &                 (maxTopChg-minTopChg)
c
                  write(11,'(a)')    '***********************'
                  write(11,'(a,3i3)') '3-Q separation scheme ', iiy, Tee, offset
                  write(11,'(a)')    '***********************'
                  call fmt_writeavg(11,actionA(:,:,:,iiy,Tee,offset))
c
                  write(12,'(a)')    '***********************'
                  write(12,'(a,3i3)') '3-Q separation scheme', iiy, Tee,offset
                  write(12,'(a)')    '***********************'
                  call fmt_writeavg(12,topChgA(:,:,:,iiy,Tee,offset))
c
               enddo               ! end iiy loop
c
            enddo                  ! end offset loop
         enddo                     ! end Tee loop
c
         close (11)
         close (12)
c
         close (9)
c
      else
c
         open (11,file=trim(reportName)//trim(fstr1)//'.correl.unf',status='unknown',form='unformatted')
         open (12,file=trim(reportName)//trim(fstr2)//'.correl.unf',status='unknown',form='unformatted')
c      
         write(11) ActionA(:,:,:,:,:,:)
         write(11) WavgA(:,:), avgActionA(:,:)
         write(12) topChgA(:,:,:,:,:,:)
         write(12) WavgA(:,:), avgtopChgA(:,:)
c      
         close(11)
         close(12)
c         
      endif
c
      end program AverageYall
