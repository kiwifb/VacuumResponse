c
c     Program to illustrate correlations between static quarks in a baryon and 
c     gluon action or topological charge
c     Deals with with T shaped Wilson loops corresponding to the Y shaped one.
c     Adapted from VacuumRespY_plan.f by F. Bissey, Oct. 2004
c
      program VacuumRespQQ

      include'VRfiles/VRCommonUse.h'
      USE L_QQLOOPS
      USE L_WRITEDQSHAPE

      IMPLICIT NONE

      include'VRfiles/VRCommonDec.h'

      character(len=132)                                         :: filename,configfile
c
c  Time averaged Quantities
c
      double precision,dimension(nx,ny,nz,nt,0:nloop,nL(that))   :: WT ! wilson loops
!HPF$ ALIGN (:,:,*,*,*,*) WITH wloops(:,:,*,*,*) :: WT

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Nomenclature for the T shape Wilson loop
c     The index running from 0 to nY correspond to the various
c     position of the three quarks in the T-loop.
c     The sets of coordinates for each value of the index are in ty-loops.f
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      double precision,dimension(0:nloop,nL(that))               :: WTavg
!HPF$ DISTRIBUTE WTavg(BLOCK,*)
      double precision,dimension(nx,ny,nz,nt)                    :: actionT
      double precision,dimension(nx,ny,nz,nt)                    :: topChgT
!HPF$ ALIGN (:,:,*,*) WITH wloops(:,:,*,*,*) :: actionT,topChgT
c
c  Correlations
c
      double precision,dimension(0:ny-1,0:nz-1,0:nt-1,0:nloop,nL(that),offmax) :: actionTC
      double precision,dimension(0:ny-1,0:nz-1,0:nt-1,0:nloop,nL(that),offmax) :: topChgTC
!HPF$ ALIGN (*,:,:,*,*,*) WITH correl(*,:,:,*,*) :: actionTC,topChgTC
c
c  Counting indexes
c
      integer                                                    :: iylp,iydir
      integer                                                    :: ix, icx, icy, icz
c
c  Begin execution
c
      write(configfile,'(a)') "./qqconfig"

      include'VRfiles/VRCommonStart.f'

      call system_clock(time1)

      call qq_loops(xhat,utr,uti,WT)

      do iylp = 1, nloop

         do it = 1, nL(that)

            WTavg(iylp,it)  = Sum( WT(:,:,:,:,iylp,it) )  / volume

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

                     do iylp = 1, nloop

                        actionTC(iy,iz,it,iylp,Tee,offset) = sum( WT(:,:,:,:,iylp,Tee) * actionTz(:,:,:,:) ) / volume

                        topChgTC(iy,iz,it,iylp,Tee,offset) = sum( WT(:,:,:,:,iylp,Tee) * topChgTz(:,:,:,:) ) / volume

                     end do     ! iylp loop end
                  end do        ! it = 0, nt-1 loop end
               end do           ! iz loop end
            end do              ! iy loop end
         end do                 ! Tee loop end
      end do                    ! offset loop end
c
      actionTC = cshift(cshift(actionTC,dim=3,shift= 4),dim=2,shift= ny/2)

      topChgTC = cshift(cshift(topChgTC,dim=3,shift= 4),dim=2,shift= ny/2)

      call system_clock(time2)
      write(*,'(a,f15.8,a)') 'time spent in the offset loop is: ',(real(time2-time1)/real(count_rate)),' s'

      filename=trim(output)//'.QQ'//trim(fstr1)
      call writeDQshape(filename,actionTC(:,:,:,:,:,:),WTavg(:,:),avgAction,offmax)

      filename=trim(output)//'.QQ'//trim(fstr2)
      call writeDQshape(filename,TopChgTC(:,:,:,:,:,:),WTavg(:,:),avgTopChg,offmax)

      include'VRfiles/VRCommonEnd.f'

      end program VacuumRespQQ

