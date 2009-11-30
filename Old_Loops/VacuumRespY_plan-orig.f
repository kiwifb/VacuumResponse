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
      character(len=132)                                         :: filename
c
c  Time averaged Quantities
c
      double precision,dimension(nx,ny,nz,nt,0:nYloop,nL(that),2)  :: WYp,WYm ! Y-shape wilson loops
!HPF$ ALIGN (:,:,*,*,*,*,*) WITH wloops(:,:,*,*,*) :: WYp,WYm

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Nomenclature for the Y shape Wilson loop
c     The index running from 0 to nY correspond to the various
c     position of the three quarks in the Y-loop.
c     The sets of coordinates for each value of the index are in y-loops.f
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      double precision,dimension(0:nYloop,nL(that),2)            :: WYpavg,WYmavg
!HPF$ DISTRIBUTE WYpavg(BLOCK,*,*)
!HPF$ DISTRIBUTE WYmavg(BLOCK,*,*)
      double precision,dimension(nx,ny,nz,nt)                    :: actionT
      double precision,dimension(nx,ny,nz,nt)                    :: topChgT
!HPF$ ALIGN (:,:,*,*) WITH wloops(:,:,*,*,*) :: actionT,topChgT
c
c  Correlations
c
      double precision,dimension(0:ny-1,0:nz-1,0:nt-1,0:nYloop,nL(that),offmax,2) :: actionYpC,actionYmC
      double precision,dimension(0:ny-1,0:nz-1,0:nt-1,0:nYloop,nL(that),offmax,2) :: topChgYpC,topChgYmc
!HPF$ ALIGN (*,:,:,*,*,*,*) WITH correl(*,:,:,*,*) :: actionYpC,topChgYpC,actionYmC,topChgYmC
c
c  Counting indexes
c
      integer                                                    :: iylp,iydir
      integer                                                    :: ix, icx, icy, icz
c
c  Begin execution
c
      include'VRfiles/VRCommonStart.f'

      call system_clock(time1)
      call y_loops_mirror(xhat,yhat,utr,uti,WYp(:,:,:,:,:,:,1),WYm(:,:,:,:,:,:,1))
      call y_loops_mirror(xhat,zhat,utr,uti,WYp(:,:,:,:,:,:,2),WYm(:,:,:,:,:,:,2))

      do iylp = 0, nYloop

         do it = 1, nL(that)
               
            WYpavg(iylp,it,1) = Sum( WYp(:,:,:,:,iylp,it,1) ) / volume
            WYmavg(iylp,it,1) = Sum( WYm(:,:,:,:,iylp,it,1) ) / volume
            WYpavg(iylp,it,2) = Sum( WYp(:,:,:,:,iylp,it,2) ) / volume
            WYmavg(iylp,it,2) = Sum( WYm(:,:,:,:,iylp,it,2) ) / volume

         end do              ! it = 2, nL(that)  loop end
      end do

      call system_clock(time2)
      write(*,'(a,f15.8,a)') 'time spent constructing Y loops is: ',(real(time2-time1)/real(count_rate)),' s'
      write(*,*) offmax
c     
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

                     do iylp = 1, nYloop

                        actionYpC(iy,iz,it,iylp,Tee,offset,1) = sum( WYp(:,:,:,:,iylp,Tee,1) * actionTz(:,:,:,:) ) / volume
                        actionYmC(iy,iz,it,iylp,Tee,offset,1) = sum( WYm(:,:,:,:,iylp,Tee,1) * actionTz(:,:,:,:) ) / volume
                        actionYpC(iy,iz,it,iylp,Tee,offset,2) = sum( WYp(:,:,:,:,iylp,Tee,2) * actionTz(:,:,:,:) ) / volume
                        actionYmC(iy,iz,it,iylp,Tee,offset,2) = sum( WYm(:,:,:,:,iylp,Tee,2) * actionTz(:,:,:,:) ) / volume

                        topChgYpC(iy,iz,it,iylp,Tee,offset,1) = sum( WYp(:,:,:,:,iylp,Tee,1) * topChgTz(:,:,:,:) ) / volume
                        topChgYmC(iy,iz,it,iylp,Tee,offset,1) = sum( WYm(:,:,:,:,iylp,Tee,1) * topChgTz(:,:,:,:) ) / volume
                        topChgYpC(iy,iz,it,iylp,Tee,offset,2) = sum( WYp(:,:,:,:,iylp,Tee,2) * topChgTz(:,:,:,:) ) / volume
                        topChgYmC(iy,iz,it,iylp,Tee,offset,2) = sum( WYm(:,:,:,:,iylp,Tee,2) * topChgTz(:,:,:,:) ) / volume

                     end do     ! iylp loop end
                  end do        ! it = 0, nt-1 loop end
               end do           ! iz loop end
            end do              ! iy loop end
         end do                 ! Tee loop end
      end do                    ! offset loop end
      call system_clock(time2)
      write(*,'(a,f15.8,a)') 'time spent in the offset loop is: ',(real(time2-time1)/real(count_rate)),' s'
      
      do ix = 0, nt-1

         icx = modulo(nt - ix , nt) 

         actionYpC(:,:,ix,1:nYloop,:,:,1) = actionYpC(:,:, ix,1:nYloop,:,:,1) + actionYmC(:,:,icx,1:nYloop,:,:,1)

         topChgYpC(:,:,ix,1:nYloop,:,:,1) = topChgYpC(:,:, ix,1:nYloop,:,:,1) + topChgYmC(:,:,icx,1:nYloop,:,:,1)

      enddo
c
      WYpavg(1:nYloop,:,1) = WYpavg(1:nYloop,:,1) + WYmavg(1:nYloop,:,1)
      
      do iy = 0, nz-1
         do iz = 0, ny-1

            icy = modulo(nz - iz , nz)
            icz = modulo(     iy , ny)

            actionYpC(iz,iy,:,1:nYloop,:,:,1) = actionYpC(iz,iy,:,1:nYloop,:,:,1) + actionYpC(icz,icy,:,1:nYloop,:,:,2)

            topChgYpC(iz,iy,:,1:nYloop,:,:,1) = topChgYpC(iz,iy,:,1:nYloop,:,:,1) + topChgYpC(icz,icy,:,1:nYloop,:,:,2)

         enddo
      enddo
c
      WYpavg(1:nYloop,:,1) = WYpavg(1:nYloop,:,1) + WYpavg(1:nYloop,:,2)
      
      do iy = 0, nz-1
         do iz = 0, ny-1
            do ix = 0, nt-1

               icx = modulo(nt - ix , nt) 
               icy = modulo(nz - iz , nz)
               icz = modulo(     iy , ny)

               actionYpC(iz,iy,ix,1:nYloop,:,:,1) = actionYpC(iz,iy,ix,1:nYloop,:,:,1) + actionYmC(icz,icy,icx,1:nYloop,:,:,2)

               topChgYpC(iz,iy,ix,1:nYloop,:,:,1) = topChgYpC(iz,iy,ix,1:nYloop,:,:,1) + topChgYmC(icz,icy,icx,1:nYloop,:,:,2)

            enddo
         enddo
      enddo
c
      WYpavg(1:nYloop,:,1) = WYpavg(1:nYloop,:,1) + WYmavg(1:nYloop,:,2)

      actionYpC(:,:,:,1:nYloop,:,:,1) = actionYpC(:,:,:,1:nYloop,:,:,1) / 4
      topChgYpC(:,:,:,1:nYloop,:,:,1) = topChgYpC(:,:,:,1:nYloop,:,:,1) / 4
      WYpavg(1:nYloop,:,1) = WYpavg(1:nYloop,:,1) / 4
      
      filename=trim(output)//'.Yavg4'//trim(fstr1)
      call writeYshape(filename,actionYpC(:,:,:,:,:,:,1),WYpavg(:,:,1),avgAction,offmax)

      filename=trim(output)//'.Yavg4'//trim(fstr2)
      call writeYshape(filename,TopChgYpC(:,:,:,:,:,:,1),WYpavg(:,:,1),avgTopChg,offmax)

      include'VRfiles/VRCommonEnd.f'

      end program VacuumRespY
      
