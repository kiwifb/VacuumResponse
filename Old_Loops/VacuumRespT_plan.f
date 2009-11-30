c
c     Program to illustrate correlations between static quarks in a baryon and 
c     gluon action or topological charge
c     Deals with with T and L shaped Wilson loops at the same time.
c     Adapted from VacuumResp.f by F. Bissey, last modified Jan. 2004
c
      program VacuumRespT

      include'VRfiles/VRCommonUse.h'
      USE L_TLOOPS
      USE L_WRITELTSHAPE
      USE L_product

      IMPLICIT NONE

      include'VRfiles/VRCommonDec.h'
      character(len=132)                                         :: configfile
      double precision,dimension(nx,ny,nz,nt,nc,nc)              :: upxr,upxi !spatial link products
      double precision,dimension(nx,ny,nz,nt,nc,nc)              :: upyr,upyi !spatial link products
      double precision,dimension(nx,ny,nz,nt,nc,nc)              :: upzr,upzi !spatial link products
!HPF$ ALIGN (*,*,:,:,*,*) WITH link(*,*,:,:,*,*) :: upxr,upxi,upyr,upyi,upzr,upzi
c
c  Time averaged Quantities
c
      double precision,dimension(nx,ny,nz,nt,nL(that),2)         :: WTp,WTm ! T-shape wilson loops
!HPF$ ALIGN (:,:,*,*,*,*) WITH wloops(:,:,*,*,*) :: WTp,WTm
      double precision,dimension(nL(that),2)                     :: WTpavg,WTmavg
!HPF$ DISTRIBUTE WTpavg(BLOCK,*)
!HPF$ DISTRIBUTE WTmavg(BLOCK,*)
c
      double precision,dimension(nx,ny,nz,nt,offmax,nL(that))    :: actionT
!HPF$ ALIGN (:,:,*,*,*,*) WITH wloops(:,:,*,*,*) :: actionT
      double precision,dimension(nx,ny,nz,nt,offmax,nL(that))    :: topChgT
!HPF$ ALIGN (:,:,*,*,*,*) WITH wloops(:,:,*,*,*) :: topChgT
c
c  Correlations
c
      double precision,dimension(0:ny-1,0:nz-1,0:nt-1,nL(that),offmax,2) :: actionTpC,actionTmC
      double precision,dimension(0:ny-1,0:nz-1,0:nt-1,nL(that),offmax,2) :: topChgTpC,topChgTmC
!HPF$ ALIGN (*,:,:,*,*,*) WITH correl(*,:,:,*,*) :: actionTpC,actionTmC,topChgTpC,topChgTmC
c
c  Counting indexes
c
      integer                                                    :: il
c
c  Begin execution
c
      write(configfile,'(a)') "/home/frbissey/SandBox/tyconfig"
      include'VRfiles/VRCommonStart.f'

c     open the files for the xy plane

      open(20,file=trim(output)//'.Tp-xy'//trim(fstr1),status='unknown',form='unformatted')
      open(21,file=trim(output)//'.Tm-xy'//trim(fstr1),status='unknown',form='unformatted')

      open(22,file=trim(output)//'.Tp-xy'//trim(fstr2),status='unknown',form='unformatted')
      open(23,file=trim(output)//'.Tm-xy'//trim(fstr2),status='unknown',form='unformatted')

c     open the files for the xz plane

      open(40,file=trim(output)//'.Tp-xz'//trim(fstr1),status='unknown',form='unformatted')
      open(41,file=trim(output)//'.Tm-xz'//trim(fstr1),status='unknown',form='unformatted')

      open(42,file=trim(output)//'.Tp-xz'//trim(fstr2),status='unknown',form='unformatted')
      open(43,file=trim(output)//'.Tm-xz'//trim(fstr2),status='unknown',form='unformatted')

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
            call t_loops(xhat,upxr,upxi,yhat,upyr,upyi,utr,uti,WTp(:,:,:,:,:,1),WTm(:,:,:,:,:,1))
            call t_loops(xhat,upxr,upxi,zhat,upzr,upzi,utr,uti,WTp(:,:,:,:,:,2),WTm(:,:,:,:,:,2))
            call system_clock(time2)
            write(*,'(a,f15.8,a)') 'time spent in shapeTLloop is: ',(Real(time2-time1)/Real(count_rate)),' s'
c
            do it = 2, nL(that)

               WTpavg(it,1) = Sum( WTp(:,:,:,:,it,1) ) / volume
               WTmavg(it,1) = Sum( WTm(:,:,:,:,it,1) ) / volume
               WTpavg(it,2) = Sum( WTp(:,:,:,:,it,2) ) / volume
               WTmavg(it,2) = Sum( WTm(:,:,:,:,it,2) ) / volume

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

                           actionTpC(iy,iz,it,Tee,offset,1) = 
     &                          sum( WTp(:,:,:,:,Tee,1) * actionTz(:,:,:,:) ) / volume
                           actionTpC(iy,iz,it,Tee,offset,2) = 
     &                          sum( WTp(:,:,:,:,Tee,2) * actionTz(:,:,:,:) ) / volume
                           actionTmC(iy,iz,it,Tee,offset,1) = 
     &                          sum( WTm(:,:,:,:,Tee,1) * actionTz(:,:,:,:) ) / volume
                           actionTmC(iy,iz,it,Tee,offset,2) = 
     &                          sum( WTm(:,:,:,:,Tee,2) * actionTz(:,:,:,:) ) / volume

                           topChgTpC(iy,iz,it,Tee,offset,1) = 
     &                          sum( WTp(:,:,:,:,Tee,1) * topChgTz(:,:,:,:) ) / volume
                           topChgTpC(iy,iz,it,Tee,offset,2) = 
     &                          sum( WTp(:,:,:,:,Tee,2) * topChgTz(:,:,:,:) ) / volume
                           topChgTmC(iy,iz,it,Tee,offset,1) = 
     &                          sum( WTm(:,:,:,:,Tee,1) * topChgTz(:,:,:,:) ) / volume
                           topChgTmC(iy,iz,it,Tee,offset,2) = 
     &                          sum( WTm(:,:,:,:,Tee,2) * topChgTz(:,:,:,:) ) / volume

                        end do  ! it = 0, nt-1 loop end
                     end do     ! iz loop end
                  end do        ! iy loop end
               end do           ! Tee loop end
            end do              ! offset loop end
            call system_clock(time2)
            write(*,'(a,f15.8,a)') 'time spent in the offset loop is: ',(real(time2-time1)/real(count_rate)),' s'
c
            call writeLTshape(20,actionTpC(:,:,:,:,:,1),WTpavg(:,1),avgAction,offmax)
            call writeLTshape(21,actionTmC(:,:,:,:,:,1),WTmavg(:,1),avgAction,offmax)

            call writeLTshape(22,topChgTpC(:,:,:,:,:,1),WTpavg(:,1),avgTopChg,offmax)
            call writeLTshape(23,topChgTmC(:,:,:,:,:,1),WTmavg(:,1),avgTopChg,offmax)

            call writeLTshape(40,actionTpC(:,:,:,:,:,2),WTpavg(:,2),avgAction,offmax)
            call writeLTshape(41,actionTmC(:,:,:,:,:,2),WTmavg(:,2),avgAction,offmax)

            call writeLTshape(42,topChgTpC(:,:,:,:,:,2),WTpavg(:,2),avgTopChg,offmax)
            call writeLTshape(43,topChgTmC(:,:,:,:,:,2),WTmavg(:,2),avgTopChg,offmax)

         end do                 ! iyq loop end

      end do                    ! ixq loop end

      do il = 20, 23

         close(il)
         close(il + 20)

      end do

      include'VRfiles/VRCommonEnd.f'

      end program VacuumRespT
