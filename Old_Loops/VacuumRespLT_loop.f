c
c     Program to illustrate correlations between static quarks in a baryon and 
c     gluon action or topological charge
c     Deals with with T and L shaped Wilson loops at the same time.
c     Adapted from VacuumResp.f by F. Bissey, last modified Jan. 2004
c
      program VacuumRespLT

      include'VRfiles/VRCommonUse.h'
      USE L_LTLOOPS
      USE L_WRITELTSHAPE
      USE L_product

      IMPLICIT NONE

      include'VRfiles/VRCommonDec.h'
      double precision,dimension(nx,ny,nz,nt,nc,nc)              :: upxr,upxi !spatial link products
      double precision,dimension(nx,ny,nz,nt,nc,nc)              :: upyr,upyi !spatial link products
!HPF$ ALIGN (*,*,:,:,*,*) WITH link(*,*,:,:,*,*) :: upxr,upxi,upyr,upyi
      double precision,dimension(nx,ny,nz,nt)                    :: actionT
      double precision,dimension(nx,ny,nz,nt)                    :: topChgT
!HPF$ ALIGN (:,:,*,*) WITH wloops(:,:,*,*,*) :: actionT,topChgT
c
c  Time averaged Quantities
c
      double precision,dimension(nx,ny,nz,nt,nL(that),4)         :: WL ! L-shape wilson loops
!HPF$ ALIGN (:,:,*,*,*,*) WITH wloops(:,:,*,*,*) :: WL
      double precision,dimension(nx,ny,nz,nt,nL(that))           :: WTp,WTm ! T-shape wilson loops
!HPF$ ALIGN (:,:,*,*,*) WITH wloops(:,:,*,*,*) :: WTp,WTm
      double precision,dimension(nL(that))                       :: WTpavg,WTmavg
!HPF$ DISTRIBUTE WTpavg(BLOCK)
!HPF$ DISTRIBUTE WTmavg(BLOCK)
      double precision,dimension(nL(that),4)                     :: WLavg
!HPF$ DISTRIBUTE WLavg(BLOCK,*)
c
c  Correlations
c
      double precision,dimension(0:ny-1,0:nz-1,0:nt-1,nL(that),offmax) :: actionTpC,actionTmC
      double precision,dimension(0:ny-1,0:nz-1,0:nt-1,nL(that),offmax) :: topChgTpC,topChgTmC
!HPF$ ALIGN (*,:,:,*,*) WITH correl(*,:,:,*,*) :: actionTpC,actionTmC,topChgTpC,topChgTmC
      double precision,dimension(0:ny-1,0:nz-1,0:nt-1,nL(that),offmax,4) :: actionLC
      double precision,dimension(0:ny-1,0:nz-1,0:nt-1,nL(that),offmax,4) :: topChgLC
!HPF$ ALIGN (*,:,:,*,*,*) WITH correl(*,:,:,*,*) :: actionLC,topChgLC
c
c    counting index
c
       integer                                                   :: il
c
c  Begin execution
c
      include'VRfiles/VRCommonStart.f'

      open(11,file=trim(output)//'.Tp'//trim(fstr1),status='unknown',form='unformatted')
      open(12,file=trim(output)//'.Tm'//trim(fstr1),status='unknown',form='unformatted')

      open(15,file=trim(output)//'.Tp'//trim(fstr2),status='unknown',form='unformatted')
      open(16,file=trim(output)//'.Tm'//trim(fstr2),status='unknown',form='unformatted')

      open(19,file=trim(output)//'.L1'//trim(fstr1),status='unknown',form='unformatted')
      open(20,file=trim(output)//'.L2'//trim(fstr1),status='unknown',form='unformatted')
      open(21,file=trim(output)//'.L3'//trim(fstr1),status='unknown',form='unformatted')
      open(22,file=trim(output)//'.L4'//trim(fstr1),status='unknown',form='unformatted')

      open(23,file=trim(output)//'.L1'//trim(fstr2),status='unknown',form='unformatted')
      open(24,file=trim(output)//'.L2'//trim(fstr2),status='unknown',form='unformatted')
      open(25,file=trim(output)//'.L3'//trim(fstr2),status='unknown',form='unformatted')
      open(26,file=trim(output)//'.L4'//trim(fstr2),status='unknown',form='unformatted')

      do ixq = 0, nL(xhat)

         call system_clock(time1)
         call product(utr,uti,xhat,upxr,upxi,ixq)
         call system_clock(time2)
         write(*,'(a,i2.2,a,f15.8,a)') 'time spent in product for ixq=',ixq,' is: ',
     &        (Real(time2-time1)/Real(count_rate)),' s'
c'
         do iyq = 0, nL(yhat)

            call product(utr,uti,yhat,upyr,upyi,iyq)
            call system_clock(time1)
            call lt_loops(xhat,upxr,upxi,yhat,upyr,upyi,utr,uti,WTp,WTm,WL)
            call system_clock(time2)
            write(*,'(a,f15.8,a)') 'time spent in shapeTLloop is: ',(Real(time2-time1)/Real(count_rate)),' s'

            do it = 2, nL(that)
               
               WTpavg(it) = Sum( WTp(:,:,:,:,it) ) / volume
               WTmavg(it) = Sum( WTm(:,:,:,:,it) ) / volume
               
               forall (il = 1:4) WLavg(it,il) = Sum( WL(:,:,:,:,it,il) ) / volume
               
            end do              ! it = 2, nL(that)  loop end

            call system_clock(time1)
            do offset = 1,offmax
               do Tee = 2*offset,nL(that)
                  do it = 1,nx

                     actionT(it,:,:,:) = 0.0d0
                     topChgT(it,:,:,:) = 0.0d0

                     do is = it+offset , it+Tee-offset

                        actionT(it,:,:,:) = actionT(it,:,:,:) + reconaction(( is - ((is - 1)/nx)*nx ),:,:,:)
                        
                        topChgT(it,:,:,:) = topChgT(it,:,:,:) + abs( topQDensity(( is - ((is - 1)/nx)*nx ),:,:,:) )
                        
                     end do     ! is loop end 
                     
                     actionT(it,:,:,:) = actionT(it,:,:,:)/(Tee - 2*offset + 1)
                     topChgT(it,:,:,:) = topChgT(it,:,:,:)/(Tee - 2*offset + 1)

                  end do        ! it = 1, nx loop end

                  avgAction(Tee,offset) = Sum( actionT(:,:,:,:) ) / volume
                  avgTopChg(Tee,offset) = Sum( topChgT(:,:,:,:) ) / volume

                  do iy = 0, ny-1

c     write(*,'(a,i2,a,i2)') 'Step ',iy,' of ',ny-1

                     if ( iy == 0 ) then

                        actionTx(:,:,:,:) = actionT(:,:,:,:)
                        topChgTx(:,:,:,:) = topChgT(:,:,:,:)

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

                           actionTpC(iy,iz,it,Tee,offset) = sum( WTp(:,:,:,:,Tee) * actionTz(:,:,:,:) ) / volume
                           actionTmC(iy,iz,it,Tee,offset) = sum( WTm(:,:,:,:,Tee) * actionTz(:,:,:,:) ) / volume

                           topChgTpC(iy,iz,it,Tee,offset) = sum( WTp(:,:,:,:,Tee) * topChgTz(:,:,:,:) ) / volume
                           topChgTmC(iy,iz,it,Tee,offset) = sum( WTm(:,:,:,:,Tee) * topChgTz(:,:,:,:) ) / volume

                           forall (il = 1:4)

                             actionLC(iy,iz,it,Tee,offset,il) = sum( WL(:,:,:,:,Tee,il) * actionTz(:,:,:,:) ) / volume

                             topChgLC(iy,iz,it,Tee,offset,il) = sum( WL(:,:,:,:,Tee,il) * topChgTz(:,:,:,:) ) / volume

                           end forall

                        end do  ! it = 0, nt-1 loop end
                     end do     ! iz loop end
                  end do        ! iy loop end
               end do           ! Tee loop end
            end do              ! offset loop end
            call system_clock(time2)
            write(*,'(a,f15.8,a)') 'time spent in the offset loop is: ',(real(time2-time1)/real(count_rate)),' s'
c'
            call writeLTshape(11,actionTpC,WTpavg,avgAction,offmax)
            call writeLTshape(12,actionTmC,WTmavg,avgAction,offmax)

            call writeLTshape(15,topChgTpC,WTpavg,avgTopChg,offmax)
            call writeLTshape(16,topChgTmC,WTmavg,avgTopChg,offmax)

            do il = 1, 4 

               call writeLTshape(18+il,actionLC(:,:,:,:,:,il),WLavg(:,il),avgAction,offmax)
               call writeLTshape(22+il,topChgLC(:,:,:,:,:,il),WLavg(:,il),avgTopChg,offmax)

            end do              ! il loop end
               
            
         end do                 ! iyq loop end
      end do                    ! ixq loop end
            
      close(11)
      close(12)
         
      close(15)
      close(16)
            
      close(19)
      close(20)
      close(21)
      close(22)

      close(23)
      close(24)
      close(25)
      close(26)

      include'VRfiles/VRCommonEnd.f'
         
      end program VacuumRespLT

      

