c
c     Program to illustrate correlations between static quarks in a baryon and 
c     gluon action or topological charge
c     Deals with with T shaped Wilson loops corresponding to the Y shaped one.
c     Adapted from VacuumRespY_plan.f by F. Bissey, Oct. 2004
c
      program VacuumRespTY

      include'VRfiles/VRCommonUse.h'
      USE L_TYLOOPS
      USE L_WRITEYSHAPE

      IMPLICIT NONE

      include'VRfiles/VRCommonDec.h'
      character(len=132)                                         :: filename
c
c  Time averaged Quantities
c
      double precision,dimension(nx,ny,nz,nt,0:nYloop,nL(that),2)  :: WTp,WTm ! T-shape wilson loops
!HPF$ ALIGN (:,:,*,*,*,*,*) WITH wloops(:,:,*,*,*) :: WTp,WTm

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Nomenclature for the T shape Wilson loop
c     The index running from 0 to nY correspond to the various
c     position of the three quarks in the T-loop.
c     The sets of coordinates for each value of the index are in ty-loops.f
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      double precision,dimension(0:nYloop,nL(that),2)            :: WTpavg,WTmavg
!HPF$ DISTRIBUTE WTpavg(BLOCK,*,*)
!HPF$ DISTRIBUTE WTmavg(BLOCK,*,*)
      double precision,dimension(nx,ny,nz,nt)                    :: actionT
      double precision,dimension(nx,ny,nz,nt)                    :: topChgT
!HPF$ ALIGN (:,:,*,*) WITH wloops(:,:,*,*,*) :: actionT,topChgT
c
c  Correlations
c
      double precision,dimension(0:ny-1,0:nz-1,0:nt-1,0:nYloop,nL(that),offmax,2) :: actionTpC,actionTmC
      double precision,dimension(0:ny-1,0:nz-1,0:nt-1,0:nYloop,nL(that),offmax,2) :: topChgTpC,topChgTmc
!HPF$ ALIGN (*,:,:,*,*,*,*) WITH correl(*,:,:,*,*) :: actionTpC,topChgTpC,actionTmC,topChgTmC
c
c  Counting indexes
c
      integer                                                    :: iylp,iydir
      integer                                                    :: ix, icx, icy, icz
      integer,dimension(nYloop)                                  :: shx
!HPF$ DISTRIBUTE shx(*)
c
c  Begin execution
c
      shx(1)=1
      shx(2)=1
      shx(3)=1
      shx(4)=2
      shx(5)=3
      shx(6)=4
      shx(7)=4
      if( nYloop .gt. 7) then
         shx(8)=4
         shx(9)=5
      endif      

      include'VRfiles/VRCommonStart.f'

      call system_clock(time1)
      call ty_loops_mirror(xhat,yhat,utr,uti,WTp(:,:,:,:,:,:,1),WTm(:,:,:,:,:,:,1))
      call ty_loops_mirror(xhat,zhat,utr,uti,WTp(:,:,:,:,:,:,2),WTm(:,:,:,:,:,:,2))

      do iylp = 0, nYloop

         do it = 1, nL(that)
               
            WTpavg(iylp,it,1) = Sum( WTp(:,:,:,:,iylp,it,1) ) / volume
            WTmavg(iylp,it,1) = Sum( WTm(:,:,:,:,iylp,it,1) ) / volume
            WTpavg(iylp,it,2) = Sum( WTp(:,:,:,:,iylp,it,2) ) / volume
            WTmavg(iylp,it,2) = Sum( WTm(:,:,:,:,iylp,it,2) ) / volume

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

                        actionTpC(iy,iz,it,iylp,Tee,offset,1) = sum( WTp(:,:,:,:,iylp,Tee,1) * actionTz(:,:,:,:) ) / volume
                        actionTmC(iy,iz,it,iylp,Tee,offset,1) = sum( WTm(:,:,:,:,iylp,Tee,1) * actionTz(:,:,:,:) ) / volume
                        actionTpC(iy,iz,it,iylp,Tee,offset,2) = sum( WTp(:,:,:,:,iylp,Tee,2) * actionTz(:,:,:,:) ) / volume
                        actionTmC(iy,iz,it,iylp,Tee,offset,2) = sum( WTm(:,:,:,:,iylp,Tee,2) * actionTz(:,:,:,:) ) / volume

                        topChgTpC(iy,iz,it,iylp,Tee,offset,1) = sum( WTp(:,:,:,:,iylp,Tee,1) * topChgTz(:,:,:,:) ) / volume
                        topChgTmC(iy,iz,it,iylp,Tee,offset,1) = sum( WTm(:,:,:,:,iylp,Tee,1) * topChgTz(:,:,:,:) ) / volume
                        topChgTpC(iy,iz,it,iylp,Tee,offset,2) = sum( WTp(:,:,:,:,iylp,Tee,2) * topChgTz(:,:,:,:) ) / volume
                        topChgTmC(iy,iz,it,iylp,Tee,offset,2) = sum( WTm(:,:,:,:,iylp,Tee,2) * topChgTz(:,:,:,:) ) / volume

                     end do     ! iylp loop end
                  end do        ! it = 0, nt-1 loop end
               end do           ! iz loop end
            end do              ! iy loop end
         end do                 ! Tee loop end
      end do                    ! offset loop end

      do ix = 0, nt-1

         icx = modulo(nt - ix , nt) 

         actionTpC(:,:,ix,1:nYloop,:,:,1) = actionTpC(:,:, ix,1:nYloop,:,:,1) + actionTmC(:,:,icx,1:nYloop,:,:,1)

         topChgTpC(:,:,ix,1:nYloop,:,:,1) = topChgTpC(:,:, ix,1:nYloop,:,:,1) + topChgTmC(:,:,icx,1:nYloop,:,:,1)

      enddo
c
      WTpavg(1:nYloop,:,1) = WTpavg(1:nYloop,:,1) + WTmavg(1:nYloop,:,1)
      
      do iy = 0, nz-1
         do iz = 0, ny-1

            icy = modulo(nz - iz , nz)
            icz = modulo(     iy , ny)

            actionTpC(iz,iy,:,1:nYloop,:,:,1) = actionTpC(iz,iy,:,1:nYloop,:,:,1) + actionTpC(icz,icy,:,1:nYloop,:,:,2)

            topChgTpC(iz,iy,:,1:nYloop,:,:,1) = topChgTpC(iz,iy,:,1:nYloop,:,:,1) + topChgTpC(icz,icy,:,1:nYloop,:,:,2)

         enddo
      enddo
c
      WTpavg(1:nYloop,:,1) = WTpavg(1:nYloop,:,1) + WTpavg(1:nYloop,:,2)
      
      do iy = 0, nz-1
         do iz = 0, ny-1
            do ix = 0, nt-1

               icx = modulo(nt - ix , nt) 
               icy = modulo(nz - iz , nz)
               icz = modulo(     iy , ny)

               actionTpC(iz,iy,ix,1:nYloop,:,:,1) = actionTpC(iz,iy,ix,1:nYloop,:,:,1) + actionTmC(icz,icy,icx,1:nYloop,:,:,2)

               topChgTpC(iz,iy,ix,1:nYloop,:,:,1) = topChgTpC(iz,iy,ix,1:nYloop,:,:,1) + topChgTmC(icz,icy,icx,1:nYloop,:,:,2)

            enddo
         enddo
      enddo
c
      WTpavg(1:nYloop,:,1) = WTpavg(1:nYloop,:,1) + WTmavg(1:nYloop,:,2)

      actionTpC(:,:,:,1:nYloop,:,:,1) = actionTpC(:,:,:,1:nYloop,:,:,1) / 4
      topChgTpC(:,:,:,1:nYloop,:,:,1) = topChgTpC(:,:,:,1:nYloop,:,:,1) / 4
      WTpavg(1:nYloop,:,1) = WTpavg(1:nYloop,:,1) / 4
      
      do iylp = 1,nYloop
         
         actionTpC(:,:,:,iylp,:,:,:) = cshift(actionTpC(:,:,:,iylp,:,:,:),dim=3,shift=shx(iylp))
         actionTmC(:,:,:,iylp,:,:,:) = cshift(actionTmC(:,:,:,iylp,:,:,:),dim=3,shift=shx(iylp))

         topChgTpC(:,:,:,iylp,:,:,:) = cshift(topChgTpC(:,:,:,iylp,:,:,:),dim=3,shift=shx(iylp))
         topChgTmC(:,:,:,iylp,:,:,:) = cshift(topChgTmC(:,:,:,iylp,:,:,:),dim=3,shift=shx(iylp)) 

      end do
         
      call system_clock(time2)
      write(*,'(a,f15.8,a)') 'time spent in the offset loop is: ',(real(time2-time1)/real(count_rate)),' s'
      filename=trim(output)//'.TYavg4'//trim(fstr1)
      call writeYshape(filename,actionTpC(:,:,:,:,:,:,1),WTpavg(:,:,1),avgAction,offmax)

      filename=trim(output)//'.TYavg4'//trim(fstr2)
      call writeYshape(filename,TopChgTpC(:,:,:,:,:,:,1),WTpavg(:,:,1),avgTopChg,offmax)

      include'VRfiles/VRCommonEnd.f'

      end program VacuumRespTY
      
