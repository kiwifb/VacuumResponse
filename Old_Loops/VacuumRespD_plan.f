c
c     Program to illustrate correlations between static quarks in a baryon and 
c     gluon action or topological charge
c     Deals with with Delta shaped loops corresponding to the Y shaped one.
c     Adapted from VacuumRespTY_plan.f by F. Bissey, Jan 2006
c
      program VacuumRespDELTA

      include'VRfiles/VRCommonUse.h'
      USE L_DELTALOOPS
      USE L_WRITEDSHAPE

      IMPLICIT NONE

      include'VRfiles/VRCommonDec.h'
      character(len=132)                                         :: filename,configfile
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
c      double precision,dimension(nx,ny,nz,nt)                    :: topChgT
!HPF$ ALIGN (:,:,*,*) WITH wloops(:,:,*,*,*) :: actionT
c,topChgT
c
c  Correlations
c
      double precision,dimension(0:ny-1,0:nz-1,0:nt-1,0:nYloop,nL(that),offmax,2) :: actionTpC,actionTmC
c      double precision,dimension(0:ny-1,0:nz-1,0:nt-1,0:nYloop,nL(that),offmax,2) :: topChgTpC,topChgTmc
!HPF$ ALIGN (*,:,:,*,*,*,*) WITH correl(*,:,:,*,*) :: actionTpC,actionTmC
c,topChgTpC,topChgTmC
c
c  Counting indexes
c
      integer                                                    :: iylp,iydir
      integer                                                    :: ix, icx, icy, icz
c      integer,dimension(nYloop)                                  :: shx
c!HPF$ DISTRIBUTE shx(*)
c
c  Begin execution
c
c      shx(1)=1
c      shx(2)=1
c      shx(3)=1
c      shx(4)=2
c      shx(5)=3
c      shx(6)=4
c      shx(7)=4
c      if( nYloop .gt. 7) then
c         shx(8)=4
c         shx(9)=5
c      endif

      write(configfile,'(a)') "/home/frbissey/SandBox/Dconfig"

      include'VRfiles/VRCommonStart.f'

      call system_clock(time1)
      call delta_loops(xhat,yhat,utr,uti,WTp(:,:,:,:,:,:,1))
c      call ty_loops_mirror(xhat,yhat,utr,uti,WTp(:,:,:,:,:,:,1),WTm(:,:,:,:,:,:,1))
c      call ty_loops_mirror(xhat,zhat,utr,uti,WTp(:,:,:,:,:,:,2),WTm(:,:,:,:,:,:,2))

      do iylp = 0, nYloop

         do it = 1, nL(that)

            WTpavg(iylp,it,1) = Sum( WTp(:,:,:,:,iylp,it,1) ) / volume
c            WTmavg(iylp,it,1) = Sum( WTm(:,:,:,:,iylp,it,1) ) / volume
c            WTpavg(iylp,it,2) = Sum( WTp(:,:,:,:,iylp,it,2) ) / volume
c            WTmavg(iylp,it,2) = Sum( WTm(:,:,:,:,iylp,it,2) ) / volume

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
c               topChgT(it,:,:,:) = 0.0d0

               do is = it+offset , it+Tee-offset

                  actionT(it,:,:,:) = actionT(it,:,:,:) + reconaction(( is - ((is - 1)/nx)*nx ),:,:,:)

c                  topChgT(it,:,:,:) = topChgT(it,:,:,:) + abs( topQDensity(( is - ((is - 1)/nx)*nx ),:,:,:) )

               end do           ! is loop end 

               actionT(it,:,:,:) = actionT(it,:,:,:)/(Tee - 2*offset + 1)
c               topChgT(it,:,:,:) = topChgT(it,:,:,:)/(Tee - 2*offset + 1)

            end do              ! it = 1, nx loop end

            avgAction(Tee,offset) = Sum( actionT(:,:,:,:) ) / volume
c            avgTopChg(Tee,offset) = Sum( topChgT(:,:,:,:) ) / volume

            do iy = 0, ny-1

c     write(*,'(a,i2,a,i2)') 'Step ',iy,' of ',ny-1

               actionTx(:,:,:,:) = cshift(actionT(:,:,:,:), dim = 2, shift = iy)
c               topChgTx(:,:,:,:) = cshift(topChgT(:,:,:,:), dim = 2, shift = iy)

               do iz = 0, nz-1

                  actionTy(:,:,:,:) = cshift(actionTx(:,:,:,:), dim = 3, shift = iz)
c                  topChgTy(:,:,:,:) = cshift(topChgTx(:,:,:,:), dim = 3, shift = iz)

                  do it = 0, nt-1

                     actionTz(:,:,:,:) = cshift(actionTy(:,:,:,:), dim = 4, shift = it)
c                     topChgTz(:,:,:,:) = cshift(topChgTy(:,:,:,:), dim = 4, shift = it)

                     do iylp = 1, nYloop

                        actionTpC(iy,iz,it,iylp,Tee,offset,1) = sum( WTp(:,:,:,:,iylp,Tee,1) * actionTz(:,:,:,:) ) / volume
c                        actionTmC(iy,iz,it,iylp,Tee,offset,1) = sum( WTm(:,:,:,:,iylp,Tee,1) * actionTz(:,:,:,:) ) / volume
c                        actionTpC(iy,iz,it,iylp,Tee,offset,2) = sum( WTp(:,:,:,:,iylp,Tee,2) * actionTz(:,:,:,:) ) / volume
c                        actionTmC(iy,iz,it,iylp,Tee,offset,2) = sum( WTm(:,:,:,:,iylp,Tee,2) * actionTz(:,:,:,:) ) / volume

c                        topChgTpC(iy,iz,it,iylp,Tee,offset,1) = sum( WTp(:,:,:,:,iylp,Tee,1) * topChgTz(:,:,:,:) ) / volume
c                        topChgTmC(iy,iz,it,iylp,Tee,offset,1) = sum( WTm(:,:,:,:,iylp,Tee,1) * topChgTz(:,:,:,:) ) / volume
c                        topChgTpC(iy,iz,it,iylp,Tee,offset,2) = sum( WTp(:,:,:,:,iylp,Tee,2) * topChgTz(:,:,:,:) ) / volume
c                        topChgTmC(iy,iz,it,iylp,Tee,offset,2) = sum( WTm(:,:,:,:,iylp,Tee,2) * topChgTz(:,:,:,:) ) / volume

                     end do     ! iylp loop end
                  end do        ! it = 0, nt-1 loop end
               end do           ! iz loop end
            end do              ! iy loop end
         end do                 ! Tee loop end
      end do                    ! offset loop end

c     do ix = 0, nt-1
c
c         icx = modulo(nt - ix , nt) 
c
c         actionTpC(:,:,ix,1:nYloop,:,:,1) = actionTpC(:,:, ix,1:nYloop,:,:,1) + actionTmC(:,:,icx,1:nYloop,:,:,1)
c
c         topChgTpC(:,:,ix,1:nYloop,:,:,1) = topChgTpC(:,:, ix,1:nYloop,:,:,1) + topChgTmC(:,:,icx,1:nYloop,:,:,1)
c
c      enddo
c
c      WTpavg(1:nYloop,:,1) = WTpavg(1:nYloop,:,1) + WTmavg(1:nYloop,:,1)

c      do iy = 0, nz-1
c         do iz = 0, ny-1

c            icy = modulo(nz - iz , nz)
c            icz = modulo(     iy , ny)

c            actionTpC(iz,iy,:,1:nYloop,:,:,1) = actionTpC(iz,iy,:,1:nYloop,:,:,1) + actionTpC(icz,icy,:,1:nYloop,:,:,2)

c            topChgTpC(iz,iy,:,1:nYloop,:,:,1) = topChgTpC(iz,iy,:,1:nYloop,:,:,1) + topChgTpC(icz,icy,:,1:nYloop,:,:,2)

c         enddo
c      enddo
c
c      WTpavg(1:nYloop,:,1) = WTpavg(1:nYloop,:,1) + WTpavg(1:nYloop,:,2)

c      do iy = 0, nz-1
c         do iz = 0, ny-1
c            do ix = 0, nt-1

c               icx = modulo(nt - ix , nt) 
c               icy = modulo(nz - iz , nz)
c               icz = modulo(     iy , ny)

c               actionTpC(iz,iy,ix,1:nYloop,:,:,1) = actionTpC(iz,iy,ix,1:nYloop,:,:,1) + actionTmC(icz,icy,icx,1:nYloop,:,:,2)

c               topChgTpC(iz,iy,ix,1:nYloop,:,:,1) = topChgTpC(iz,iy,ix,1:nYloop,:,:,1) + topChgTmC(icz,icy,icx,1:nYloop,:,:,2)

c            enddo
c         enddo
c      enddo
c
c      WTpavg(1:nYloop,:,1) = WTpavg(1:nYloop,:,1) + WTmavg(1:nYloop,:,2)

c      actionTpC(:,:,:,1:nYloop,:,:,1) = actionTpC(:,:,:,1:nYloop,:,:,1) / 4
c      topChgTpC(:,:,:,1:nYloop,:,:,1) = topChgTpC(:,:,:,1:nYloop,:,:,1) / 4
c      WTpavg(1:nYloop,:,1) = WTpavg(1:nYloop,:,1) / 4

c      do iylp = 1,nYloop

c         actionTpC(:,:,:,iylp,:,:,:) = cshift(actionTpC(:,:,:,iylp,:,:,:),dim=3,shift=shx(iylp))
c         actionTmC(:,:,:,iylp,:,:,:) = cshift(actionTmC(:,:,:,iylp,:,:,:),dim=3,shift=shx(iylp))

c         topChgTpC(:,:,:,iylp,:,:,:) = cshift(topChgTpC(:,:,:,iylp,:,:,:),dim=3,shift=shx(iylp))
c         topChgTmC(:,:,:,iylp,:,:,:) = cshift(topChgTmC(:,:,:,iylp,:,:,:),dim=3,shift=shx(iylp)) 

c      end do

      call system_clock(time2)
      write(*,'(a,f15.8,a)') 'time spent in the offset loop is: ',(real(time2-time1)/real(count_rate)),' s'
      filename=trim(output)//'.Davg4'//trim(fstr1)

      open (11,file=trim(filename),status='unknown')

      do Tee = 2,nL(that)
         do offset = 1,offmax
c
            if ( Tee - 2*offset .lt. 0 ) then
               cycle
            endif
c
            do iylp = 1, nYloop
c
               write(11,'(a)')    '*********************'
               write(11,'(a,3i3)') 'Delta Loop ', iylp, Tee,offset
               write(11,'(a)')    '*********************'
               write(11,'(f15.8)') actionTpC(:,:,:,iylp,Tee,offset,1) / (WTpavg(iylp,Tee,1) *avgAction(Tee,offset))
c
               enddo               ! end iylp loop
c
            enddo                  ! end offset loop
         enddo                     ! end Tee loop
c
         close (11)

c      call writeDshape(filename,actionTpC(:,:,:,:,:,:,1),WTpavg(:,:,1),avgAction,offmax)

c      filename=trim(output)//'.TYavg4'//trim(fstr2)
c      call writeYshape(filename,TopChgTpC(:,:,:,:,:,:,1),WTpavg(:,:,1),avgTopChg,offmax)

      include'VRfiles/VRCommonEnd.f'

      end program VacuumRespDelta