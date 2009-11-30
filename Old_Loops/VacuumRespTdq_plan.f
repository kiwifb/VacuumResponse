c
c     Program to illustrate correlations between static quarks in a baryon and 
c     gluon action or topological charge
c     Deals with with T shaped Wilson loops corresponding to the Y shaped one.
c     Adapted from VacuumRespY_plan.f by F. Bissey, Oct. 2004
c
      program VacuumRespTdq

      include'VRfiles/VRCommonUse.h'
      USE L_DQLOOPS
      USE L_WRITEDQSHAPE

      IMPLICIT NONE

      include'VRfiles/VRCommonDec.h'

      character(len=132)                                         :: filename,configfile
c
c  Time averaged Quantities
c
      double precision,dimension(nx,ny,nz,nt,0:ldq,nL(that),2)   :: WTp,WTm ! T-shape wilson loops
!HPF$ ALIGN (:,:,*,*,*,*,*) WITH wloops(:,:,*,*,*) :: WTp,WTm
      double precision,dimension(nx,ny,nz,nt,0:ldq,nL(that),2)   :: WT4p,WT4m ! T-shape wilson loops
!HPF$ ALIGN (:,:,*,*,*,*,*) WITH wloops(:,:,*,*,*) :: WT4p,WT4m

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Nomenclature for the T shape Wilson loop
c     The index running from 0 to nY correspond to the various
c     position of the three quarks in the T-loop.
c     The sets of coordinates for each value of the index are in ty-loops.f
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      double precision,dimension(0:ldq,nL(that),2)               :: WTpavg,WTmavg
!HPF$ DISTRIBUTE WTpavg(BLOCK,*,*)
!HPF$ DISTRIBUTE WTmavg(BLOCK,*,*)
      double precision,dimension(0:ldq,nL(that),2)               :: WT4pavg,WT4mavg
!HPF$ DISTRIBUTE WT4pavg(BLOCK,*,*)
!HPF$ DISTRIBUTE WT4mavg(BLOCK,*,*)
      double precision,dimension(nx,ny,nz,nt)                    :: actionT
      double precision,dimension(nx,ny,nz,nt)                    :: topChgT
!HPF$ ALIGN (:,:,*,*) WITH wloops(:,:,*,*,*) :: actionT,topChgT
c
c  Correlations
c
      double precision,dimension(0:ny-1,0:nz-1,0:nt-1,0:ldq,nL(that),offmax,2) :: actionTpC,actionTmC
      double precision,dimension(0:ny-1,0:nz-1,0:nt-1,0:ldq,nL(that),offmax,2) :: topChgTpC,topChgTmc
!HPF$ ALIGN (*,:,:,*,*,*,*) WITH correl(*,:,:,*,*) :: actionTpC,topChgTpC,actionTmC,topChgTmC
      double precision,dimension(0:ny-1,0:nz-1,0:nt-1,0:ldq,nL(that),offmax,2) :: actionT4pC,actionT4mC
      double precision,dimension(0:ny-1,0:nz-1,0:nt-1,0:ldq,nL(that),offmax,2) :: topChgT4pC,topChgT4mc
!HPF$ ALIGN (*,:,:,*,*,*,*) WITH correl(*,:,:,*,*) :: actionT4pC,topChgT4pC,actionT4mC,topChgT4mC
c
c  Counting indexes
c
      integer                                                    :: iylp,iydir
      integer                                                    :: ix, icx, icy, icz
c
c  Begin execution
c
      write(configfile,'(a)') "./dqconfig"

      include'VRfiles/VRCommonStart.f'

      call system_clock(time1)
      call t_loops(xhat,yhat,utr,uti,
     &             WTp(:,:,:,:,:,:,1),WTm(:,:,:,:,:,:,1),WT4p(:,:,:,:,:,:,1),WT4m(:,:,:,:,:,:,1))
      call t_loops(xhat,zhat,utr,uti,
     &             WTp(:,:,:,:,:,:,2),WTm(:,:,:,:,:,:,2),WT4p(:,:,:,:,:,:,2),WT4m(:,:,:,:,:,:,2))

      do iylp = 1, ldq

         do it = 1, nL(that)

            WTpavg(iylp,it,1)  = Sum( WTp(:,:,:,:,iylp,it,1) )  / volume
            WTmavg(iylp,it,1)  = Sum( WTm(:,:,:,:,iylp,it,1) )  / volume
            WTpavg(iylp,it,2)  = Sum( WTp(:,:,:,:,iylp,it,2) )  / volume
            WTmavg(iylp,it,2)  = Sum( WTm(:,:,:,:,iylp,it,2) )  / volume
            WT4pavg(iylp,it,1) = Sum( WT4p(:,:,:,:,iylp,it,1) ) / volume
            WT4mavg(iylp,it,1) = Sum( WT4m(:,:,:,:,iylp,it,1) ) / volume
            WT4pavg(iylp,it,2) = Sum( WT4p(:,:,:,:,iylp,it,2) ) / volume
            WT4mavg(iylp,it,2) = Sum( WT4m(:,:,:,:,iylp,it,2) ) / volume

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

                     do iylp = 1, ldq

                        actionTpC(iy,iz,it,iylp,Tee,offset,1) = sum( WTp(:,:,:,:,iylp,Tee,1) * actionTz(:,:,:,:) ) / volume
                        actionTmC(iy,iz,it,iylp,Tee,offset,1) = sum( WTm(:,:,:,:,iylp,Tee,1) * actionTz(:,:,:,:) ) / volume
                        actionTpC(iy,iz,it,iylp,Tee,offset,2) = sum( WTp(:,:,:,:,iylp,Tee,2) * actionTz(:,:,:,:) ) / volume
                        actionTmC(iy,iz,it,iylp,Tee,offset,2) = sum( WTm(:,:,:,:,iylp,Tee,2) * actionTz(:,:,:,:) ) / volume

                        topChgTpC(iy,iz,it,iylp,Tee,offset,1) = sum( WTp(:,:,:,:,iylp,Tee,1) * topChgTz(:,:,:,:) ) / volume
                        topChgTmC(iy,iz,it,iylp,Tee,offset,1) = sum( WTm(:,:,:,:,iylp,Tee,1) * topChgTz(:,:,:,:) ) / volume
                        topChgTpC(iy,iz,it,iylp,Tee,offset,2) = sum( WTp(:,:,:,:,iylp,Tee,2) * topChgTz(:,:,:,:) ) / volume
                        topChgTmC(iy,iz,it,iylp,Tee,offset,2) = sum( WTm(:,:,:,:,iylp,Tee,2) * topChgTz(:,:,:,:) ) / volume

                        actionT4pC(iy,iz,it,iylp,Tee,offset,1) = sum( WT4p(:,:,:,:,iylp,Tee,1) * actionTz(:,:,:,:) ) / volume
                        actionT4mC(iy,iz,it,iylp,Tee,offset,1) = sum( WT4m(:,:,:,:,iylp,Tee,1) * actionTz(:,:,:,:) ) / volume
                        actionT4pC(iy,iz,it,iylp,Tee,offset,2) = sum( WT4p(:,:,:,:,iylp,Tee,2) * actionTz(:,:,:,:) ) / volume
                        actionT4mC(iy,iz,it,iylp,Tee,offset,2) = sum( WT4m(:,:,:,:,iylp,Tee,2) * actionTz(:,:,:,:) ) / volume

                        topChgT4pC(iy,iz,it,iylp,Tee,offset,1) = sum( WT4p(:,:,:,:,iylp,Tee,1) * topChgTz(:,:,:,:) ) / volume
                        topChgT4mC(iy,iz,it,iylp,Tee,offset,1) = sum( WT4m(:,:,:,:,iylp,Tee,1) * topChgTz(:,:,:,:) ) / volume
                        topChgT4pC(iy,iz,it,iylp,Tee,offset,2) = sum( WT4p(:,:,:,:,iylp,Tee,2) * topChgTz(:,:,:,:) ) / volume
                        topChgT4mC(iy,iz,it,iylp,Tee,offset,2) = sum( WT4m(:,:,:,:,iylp,Tee,2) * topChgTz(:,:,:,:) ) / volume

                     end do     ! iylp loop end
                  end do        ! it = 0, nt-1 loop end
               end do           ! iz loop end
            end do              ! iy loop end
         end do                 ! Tee loop end
      end do                    ! offset loop end

      do ix = 0, nt-1

         icx = modulo(nt - ix , nt) 

         actionTpC(:,:,ix,1:ldq,:,:,1) = actionTpC(:,:, ix,1:ldq,:,:,1) + actionTmC(:,:,icx,1:ldq,:,:,1)

         topChgTpC(:,:,ix,1:ldq,:,:,1) = topChgTpC(:,:, ix,1:ldq,:,:,1) + topChgTmC(:,:,icx,1:ldq,:,:,1)

         actionT4pC(:,:,ix,1:ldq,:,:,1) = actionT4pC(:,:, ix,1:ldq,:,:,1) + actionT4mC(:,:,icx,1:ldq,:,:,1)

         topChgT4pC(:,:,ix,1:ldq,:,:,1) = topChgT4pC(:,:, ix,1:ldq,:,:,1) + topChgT4mC(:,:,icx,1:ldq,:,:,1)

      end do
c
      WTpavg(1:ldq,:,1) = WTpavg(1:ldq,:,1) + WTmavg(1:ldq,:,1)
      WT4pavg(1:ldq,:,1) = WT4pavg(1:ldq,:,1) + WT4mavg(1:ldq,:,1)

      do iy = 0, nz-1
         do iz = 0, ny-1

            icy = modulo(nz - iz , nz)
            icz = modulo(     iy , ny)

            actionTpC(iz,iy,:,1:ldq,:,:,1) = actionTpC(iz,iy,:,1:ldq,:,:,1) + actionTpC(icz,icy,:,1:ldq,:,:,2)

            topChgTpC(iz,iy,:,1:ldq,:,:,1) = topChgTpC(iz,iy,:,1:ldq,:,:,1) + topChgTpC(icz,icy,:,1:ldq,:,:,2)

            actionT4pC(iz,iy,:,1:ldq,:,:,1) = actionT4pC(iz,iy,:,1:ldq,:,:,1) + actionT4pC(icz,icy,:,1:ldq,:,:,2)

            topChgT4pC(iz,iy,:,1:ldq,:,:,1) = topChgT4pC(iz,iy,:,1:ldq,:,:,1) + topChgT4pC(icz,icy,:,1:ldq,:,:,2)

         end do
      end do
c
      WTpavg(1:ldq,:,1) = WTpavg(1:ldq,:,1) + WTpavg(1:ldq,:,2)
      WT4pavg(1:ldq,:,1) = WT4pavg(1:ldq,:,1) + WT4pavg(1:ldq,:,2)

      do iy = 0, nz-1
         do iz = 0, ny-1
            do ix = 0, nt-1

               icx = modulo(nt - ix , nt) 
               icy = modulo(nz - iz , nz)
               icz = modulo(     iy , ny)

               actionTpC(iz,iy,ix,1:ldq,:,:,1) = actionTpC(iz,iy,ix,1:ldq,:,:,1) + actionTmC(icz,icy,icx,1:ldq,:,:,2)

               topChgTpC(iz,iy,ix,1:ldq,:,:,1) = topChgTpC(iz,iy,ix,1:ldq,:,:,1) + topChgTmC(icz,icy,icx,1:ldq,:,:,2)

               actionT4pC(iz,iy,ix,1:ldq,:,:,1) = actionT4pC(iz,iy,ix,1:ldq,:,:,1) + actionT4mC(icz,icy,icx,1:ldq,:,:,2)

               topChgT4pC(iz,iy,ix,1:ldq,:,:,1) = topChgT4pC(iz,iy,ix,1:ldq,:,:,1) + topChgT4mC(icz,icy,icx,1:ldq,:,:,2)

            end do
         end do
      end do
c
      WTpavg(1:ldq,:,1) = WTpavg(1:ldq,:,1) + WTmavg(1:ldq,:,2)
      WT4pavg(1:ldq,:,1) = WT4pavg(1:ldq,:,1) + WT4mavg(1:ldq,:,2)

      actionTpC(:,:,:,1:ldq,:,:,1) = actionTpC(:,:,:,1:ldq,:,:,1) / 4
      topChgTpC(:,:,:,1:ldq,:,:,1) = topChgTpC(:,:,:,1:ldq,:,:,1) / 4
      WTpavg(1:ldq,:,1) = WTpavg(1:ldq,:,1) / 4

      actionT4pC(:,:,:,1:ldq,:,:,1) = actionT4pC(:,:,:,1:ldq,:,:,1) / 4
      topChgT4pC(:,:,:,1:ldq,:,:,1) = topChgT4pC(:,:,:,1:ldq,:,:,1) / 4
      WT4pavg(1:ldq,:,1) = WT4pavg(1:ldq,:,1) / 4

      actionTpC = cshift(cshift(actionTpC,dim=3,shift= 4),dim=2,shift= ny/2)
      actionT4pC = cshift(cshift(actionT4pC,dim=3,shift= 4),dim=2,shift= ny/2)

      topChgTpC = cshift(cshift(topChgTpC,dim=3,shift= 4),dim=2,shift= ny/2)
      topChgT4pC = cshift(cshift(topChgT4pC,dim=3,shift= 4),dim=2,shift= ny/2)

      call system_clock(time2)
      write(*,'(a,f15.8,a)') 'time spent in the offset loop is: ',(real(time2-time1)/real(count_rate)),' s'

      filename=trim(output)//'.DQ2avg4'//trim(fstr1)
      call writeDQshape(filename,actionTpC(:,:,:,:,:,:,1),WTpavg(:,:,1),avgAction,offmax)

      filename=trim(output)//'.DQ2avg4'//trim(fstr2)
      call writeDQshape(filename,TopChgTpC(:,:,:,:,:,:,1),WTpavg(:,:,1),avgTopChg,offmax)

      filename=trim(output)//'.DQ4avg4'//trim(fstr1)
      call writeDQshape(filename,actionT4pC(:,:,:,:,:,:,1),WT4pavg(:,:,1),avgAction,offmax)

      filename=trim(output)//'.DQ4avg4'//trim(fstr2)
      call writeDQshape(filename,TopChgT4pC(:,:,:,:,:,:,1),WT4pavg(:,:,1),avgTopChg,offmax)

      include'VRfiles/VRCommonEnd.f'

      end program VacuumRespTdq

