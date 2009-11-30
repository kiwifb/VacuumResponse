c
c     Program to illustrate correlations between static quarks in a baryon and 
c     gluon action or topological charge
c     Deals with with 8 L shaped Wilson loops, 4 in each possible plan.
c     Adapted from VacuumResp.f by F. Bissey.
c     FB050324
c
      program VacuumRespL

      include'VRfiles/VRCommonUse.h'
      USE L_LLOOPS
      USE L_WRITELTSHAPE
      USE L_product

      IMPLICIT NONE

      include'VRfiles/VRCommonDec.h'
      double precision,dimension(nx,ny,nz,nt,nc,nc)              :: upxr,upxi !spatial link products
      double precision,dimension(nx,ny,nz,nt,nc,nc)              :: upyr,upyi !spatial link products
      double precision,dimension(nx,ny,nz,nt,nc,nc)              :: upzr,upzi !spatial link products
!HPF$ ALIGN (*,*,:,:,*,*) WITH link(*,*,:,:,*,*) :: upxr,upxi,upyr,upyi,upzr,upzi
c
c  Time averaged Quantities
c
      double precision,dimension(nx,ny,nz,nt,nL(that),4,2)       :: WL ! L-shape wilson loops
!HPF$ ALIGN (:,:,*,*,*,*,*) WITH wloops(:,:,*,*,*) :: WL
      double precision,dimension(nL(that),4,2)                   :: WLavg
!HPF$ DISTRIBUTE WLavg(BLOCK,*,*)
c
      double precision,dimension(nx,ny,nz,nt,offmax,nL(that))    :: actionT
!HPF$ ALIGN (:,:,*,*,*,*) WITH wloops(:,:,*,*,*) :: actionT
      double precision,dimension(nx,ny,nz,nt,offmax,nL(that))    :: topChgT
!HPF$ ALIGN (:,:,*,*,*,*) WITH wloops(:,:,*,*,*) :: topChgT
c
c  Correlations
c
      double precision,dimension(0:ny-1,0:nz-1,0:nt-1,nL(that),offmax,4,2) :: actionLC
      double precision,dimension(0:ny-1,0:nz-1,0:nt-1,nL(that),offmax,4,2) :: topChgLC
!HPF$ ALIGN (*,:,:,*,*,*,*) WITH correl(*,:,:,*,*) :: actionLC,topChgLC
c
c  Counting indexes
c
      integer                                                    :: il
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
      
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     initialize the epsilon index for the calls from lt_oop
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

         do istep = 1,10
            
            call ape_smear3D(utr,uti,alpha,that)

         end do

         if ( Correlate == 1) then
               
            fstr1='.action.correl.unf'
            fstr2='.topChg.correl.unf'
            
         else

            fstr1='.ele.correl.unf'
            fstr2='.mag.correl.unf'
            
         endif

c     open the files for the xy plane

      open(24,file=trim(output)//'.L1-xy'//trim(fstr1),status='unknown',form='unformatted')
      open(25,file=trim(output)//'.L2-xy'//trim(fstr1),status='unknown',form='unformatted')
      open(26,file=trim(output)//'.L3-xy'//trim(fstr1),status='unknown',form='unformatted')
      open(27,file=trim(output)//'.L4-xy'//trim(fstr1),status='unknown',form='unformatted')

      open(28,file=trim(output)//'.L1-xy'//trim(fstr2),status='unknown',form='unformatted')
      open(29,file=trim(output)//'.L2-xy'//trim(fstr2),status='unknown',form='unformatted')
      open(30,file=trim(output)//'.L3-xy'//trim(fstr2),status='unknown',form='unformatted')
      open(31,file=trim(output)//'.L4-xy'//trim(fstr2),status='unknown',form='unformatted')

c     open the files for the xz plane

      open(44,file=trim(output)//'.L1-xz'//trim(fstr1),status='unknown',form='unformatted')
      open(45,file=trim(output)//'.L2-xz'//trim(fstr1),status='unknown',form='unformatted')
      open(46,file=trim(output)//'.L3-xz'//trim(fstr1),status='unknown',form='unformatted')
      open(47,file=trim(output)//'.L4-xz'//trim(fstr1),status='unknown',form='unformatted')

      open(48,file=trim(output)//'.L1-xz'//trim(fstr2),status='unknown',form='unformatted')
      open(49,file=trim(output)//'.L2-xz'//trim(fstr2),status='unknown',form='unformatted')
      open(50,file=trim(output)//'.L3-xz'//trim(fstr2),status='unknown',form='unformatted')
      open(51,file=trim(output)//'.L4-xz'//trim(fstr2),status='unknown',form='unformatted')

      do offset = 1,offmax
         do Tee = 2*offset,nL(that)
            do it = 1,nx

               actionT(it,:,:,:,offset,Tee) = 0.0d0
               topChgT(it,:,:,:,offset,Tee) = 0.0d0

               do is = it+offset , it+Tee-offset

                  actionT(it,:,:,:,offset,Tee) = actionT(it,:,:,:,offset,Tee) + 
     &                 reconaction(( is - ((is - 1)/nx)*nx ),:,:,:)
                           
                  topChgT(it,:,:,:,offset,Tee) = topChgT(it,:,:,:,offset,Tee) + 
     &                 abs( topQDensity(( is - ((is - 1)/nx)*nx ),:,:,:) )
                  
               end do           ! is loop end 

               actionT(it,:,:,:,offset,Tee) = 
     &              actionT(it,:,:,:,offset,Tee)/(Tee - 2*offset + 1)
               topChgT(it,:,:,:,offset,Tee) = 
     &              topChgT(it,:,:,:,offset,Tee)/(Tee - 2*offset + 1)

            end do              ! it = 1, nx loop end

            avgAction(Tee,offset) = Sum( actionT(:,:,:,:,offset,Tee) ) / volume
            avgTopChg(Tee,offset) = Sum( topChgT(:,:,:,:,offset,Tee) ) / volume

         end do                 ! Tee loop end 
      end do                    ! offset loop end

      do ixq = 0, nL(xhat)

         call system_clock(time1)
         call product(utr,uti,xhat,upxr,upxi,ixq)
         call system_clock(time2)
         write(*,'(a,i2.2,a,f15.8,a)') 'time spent in product for ixq=',ixq,' is: ',
     &        (Real(time2-time1)/Real(count_rate)),' s'
c            
         do iyq = 0, nL(yhat)

            call product(utr,uti,yhat,upyr,upyi,iyq)
            call product(utr,uti,zhat,upzr,upzi,iyq)
            call system_clock(time1)
            call l_loops(xhat,upxr,upxi,yhat,upyr,upyi,utr,uti,WL(:,:,:,:,:,:,1))
            call l_loops(xhat,upxr,upxi,zhat,upzr,upzi,utr,uti,WL(:,:,:,:,:,:,2))
            call system_clock(time2)
            write(*,'(a,f15.8,a)') 'time spent in shapeTLloop is: ',(Real(time2-time1)/Real(count_rate)),' s'
c'
            do it = 2, nL(that)
               
               forall (il = 1:4) WLavg(it,il,1) = Sum( WL(:,:,:,:,it,il,1) ) / volume
               forall (il = 1:4) WLavg(it,il,2) = Sum( WL(:,:,:,:,it,il,2) ) / volume

            end do              ! it = 2, nL(that)  loop end

            call system_clock(time1)
            do offset = 1,offmax
               do Tee = 2*offset,nL(that)
                  do iy = 0, ny-1

c     write(*,'(a,i2,a,i2)') 'Step ',iy,' of ',ny-1

                     if ( iy == 0 ) then

                        actionTx(:,:,:,:) = actionT(:,:,:,:,offset,Tee)
                        topChgTx(:,:,:,:) = topChgT(:,:,:,:,offset,Tee)

                     else
                           
                        actionTx(:,:,:,:) = cshift(actionTx(:,:,:,:), shift = 1, dim =2)
                        topChgTx(:,:,:,:) = cshift(topChgTx(:,:,:,:), shift = 1, dim =2)

                     end if
                        
                     do iz = 0, nz-1
                           
                        if ( iz == 0 ) then
                           
                           actionTy(:,:,:,:) = actionTx(:,:,:,:)
                           topChgTy(:,:,:,:) = topChgTx(:,:,:,:)

                        else

                           actionTy(:,:,:,:) = cshift(actionTy(:,:,:,:), shift = 1, dim = 3)
                           topChgTy(:,:,:,:) = cshift(topChgTy(:,:,:,:), shift = 1, dim = 3)

                        end if

                        do it = 0, nt-1
                           
                           if ( it == 0 ) then
                                 
                              actionTz(:,:,:,:) = actionTy(:,:,:,:)
                              topChgTz(:,:,:,:) = topChgTy(:,:,:,:)

                           else
                                 
                              actionTz(:,:,:,:) = cshift(actionTz(:,:,:,:), shift = 1, dim = 4)
                              topChgTz(:,:,:,:) = cshift(topChgTz(:,:,:,:), shift = 1, dim = 4)
                
                           end if

                           forall (il = 1:4)

                             actionLC(iy,iz,it,Tee,offset,il,1) = 
     &                          sum( WL(:,:,:,:,Tee,il,1) * actionTz(:,:,:,:) ) / volume
                             actionLC(iy,iz,it,Tee,offset,il,2) = 
     &                            sum( WL(:,:,:,:,Tee,il,2) * actionTz(:,:,:,:) ) / volume
                             topChgLC(iy,iz,it,Tee,offset,il,1) = 
     &                            sum( WL(:,:,:,:,Tee,il,1) * topChgTz(:,:,:,:) ) / volume
                             topChgLC(iy,iz,it,Tee,offset,il,2) = 
     &                            sum( WL(:,:,:,:,Tee,il,2) * topChgTz(:,:,:,:) ) / volume

                           end forall

                        end do  ! it = 0, nt-1 loop end
                     end do     ! iz loop end
                  end do        ! iy loop end
               end do           ! Tee loop end
            end do              ! offset loop end
            call system_clock(time2)
            write(*,'(a,f15.8,a)') 'time spent in the offset loop is: ',(real(time2-time1)/real(count_rate)),' s'
c'
            do il = 1, 4 

               call writeLTshape(23+il,actionLC(:,:,:,:,:,il,1),WLavg(:,il,1),avgAction,offmax)
               call writeLTshape(27+il,topChgLC(:,:,:,:,:,il,1),WLavg(:,il,1),avgTopChg,offmax)

               call writeLTshape(43+il,actionLC(:,:,:,:,:,il,2),WLavg(:,il,2),avgAction,offmax)
               call writeLTshape(47+il,topChgLC(:,:,:,:,:,il,2),WLavg(:,il,2),avgTopChg,offmax)

            end do              ! il loop end

         end do                 ! iyq loop end
         
      end do                    ! ixq loop end
      
      do il = 23, 31
            
         close(il)
         close(il + 20)

      end do

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

      end program VacuumRespL

      

