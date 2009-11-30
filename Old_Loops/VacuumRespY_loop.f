c
c     Program to illustrate correlations between static quarks in a baryon and 
c     gluon action or topological charge
c     Deals with with Y shaped Wilson loops.
c     Adapted from VacuumRespLT.f by F. Bissey, Jan. 2004
c
      program VacuumRespY

      include'VRfiles/VRCommonUse.h'
      USE L_YLOOPS
      USE L_WRITEYSHAPE

      IMPLICIT NONE

      include'VRfiles/VRCommonDec.h'
      character(len=132)                                         :: filename,configfile
c
c  Time averaged Quantities
c
      double precision,dimension(nx,ny,nz,nt,0:nYloop,nL(that))  :: WY ! Y-shape wilson loops
!HPF$ ALIGN (:,:,*,*,*,*) WITH wloops(:,:,*,*,*) :: WY

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Nomenclature for the Y shape Wilson loop
c     The index running from 0 to nY correspond to the various
c     position of the three quarks in the Y-loop.
c     The sets of coordinates for each loop is in y-loops.f
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      double precision,dimension(0:nYloop,nL(that))              :: WYavg
!HPF$ DISTRIBUTE WYavg(BLOCK,*)
      double precision,dimension(nx,ny,nz,nt)                    :: actionT
      double precision,dimension(nx,ny,nz,nt)                    :: topChgT
!HPF$ ALIGN (:,:,*,*) WITH wloops(:,:,*,*,*) :: actionT,topChgT
c
c  Correlations
c
      double precision,dimension(0:ny-1,0:nz-1,0:nt-1,0:nYloop,nL(that),offmax) :: actionYC
      double precision,dimension(0:ny-1,0:nz-1,0:nt-1,0:nYloop,nL(that),offmax) :: topChgYC
!HPF$ ALIGN (*,:,:,*,*,*) WITH correl(*,:,:,*,*) :: actionYC,topChgYC
c
c  Counting indexes
c
      integer                                                    :: iiy
c
c  Begin execution
c
      write(configfile,'(a)') "/home/frbissey/SandBox/yconfig"
      include'VRfiles/VRCommonStart.f'

      call system_clock(time1)
      call y_loops(xhat,yhat,utr,uti,WY)
         
      do iiy = 0, nYloop

         do it = 2, nL(that)
               
            WYavg(iiy,it) = Sum( WY(:,:,:,:,iiy,it) ) / volume

         end do                 ! it = 2, nL(that)  loop end
      end do
      call system_clock(time2)
      write(*,'(a,f15.8,a)') 'time spent constructing Y loops is: ',(real(time2-time1)/real(count_rate)),' s'
      write(*,*) offmax

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
      filename=trim(output)//'.Y'//trim(fstr1)
      call writeYshape(filename,actionYC(:,:,:,:,:,:),WYavg(:,:),avgAction,offmax)
      
      filename=trim(output)//'.Y'//trim(fstr2)
      call writeYshape(filename,TopChgYC(:,:,:,:,:,:),WYavg(:,:),avgTopChg,offmax)
      
      include'VRfiles/VRCommonEnd.f'

      end program VacuumRespY

      

