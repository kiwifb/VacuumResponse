!
!     ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!     subroutine pseudosweep updates the current link using a psedo-heatbath
!     algorithm. It calls the staples and pseudoheat.
!
!     Author: Frederi! D.R. Bonnet & Derek B. Leinweber
!     Date: 4th of August 1998
!
!     Modified: 4th January 2002 to replace subsequent link times
!       staples multiplies by su2 subgroup multiplies.
!       Derek B. Leinweber
!
      MODULE GS_PSEUDOSWEEP

      CONTAINS

      subroutine pseudosweep(ur,ui,mask,umask,beta,itype,uzero)

      USE GS_LATTICESIZE
      USE GS_STAPLES
      USE GS_PSEUDOHEAT

      implicit none
!     global variables

      integer,parameter                                       :: nc=3,ncsu2=2
      integer,parameter                                       :: nsigma=ncsu2*ncsu2
      integer,parameter                                       :: mu=4
      integer,parameter                                       :: nmask=16

      integer                                                 :: umask
      integer                                                 :: itype

      double precision                                        :: uzero
      double precision                                        :: beta

      double precision,dimension(nx,ny,nz,nt,mu,nc,nc)        :: ur,ui
!HPF$ DISTRIBUTE ur(*,*,BLOCK,BLOCK,*,*,*)
!HPF$ DISTRIBUTE ui(*,*,BLOCK,BLOCK,*,*,*)
      logical,dimension(nx,ny,nz,nt,mu,nmask)                 :: mask
!HPF$ DISTRIBUTE mask(*,*,BLOCK,BLOCK,*,*)

!     local variables

      integer,parameter                                       :: nloop=2

      double precision,dimension(nx,ny,nz,nt,nc,nc)           :: stapler,staplei
      double precision,dimension(nx,ny,nz,nt,nc,nc)           :: ltsr,ltsi
      double precision,dimension(nx,ny,nz,nt,nc,nc)           :: urprmsu3,uiprmsu3
      double precision,dimension(nx,ny,nz,nt,nc,nc)           :: ltsrprm,ltsiprm
!HPF$ DISTRIBUTE stapler(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE staplei(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE ltsr(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE ltsi(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE urprmsu3(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE uiprmsu3(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE ltsrprm(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE ltsiprm(*,*,BLOCK,BLOCK,*,*)
      double precision,dimension(nx,ny,nz,nt,ncsu2,ncsu2)     :: phbsr,phbsi
      double precision,dimension(nx,ny,nz,nt,ncsu2,ncsu2)     :: ursu2,uisu2
!HPF$ DISTRIBUTE phbsr(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE phbsi(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE ursu2(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE uisu2(*,*,BLOCK,BLOCK,*,*)

      double precision                                        :: betanew
      integer                                                 :: imask,ihat,ic,jc,ic3
      integer                                                 :: iloop

!
!  Timer Support
!
      INTEGER start_count, end_count, count_rate
      REAL    elapsed_time
!
!  Useage
!
!  CALL SYSTEM_CLOCK(start_count, count_rate)
!  CALL SYSTEM_CLOCK(end_count)
!  elapsed_time = REAL(end_count - start_count) / REAL(count_rate)
!

!     starts of the running commands

      CALL SYSTEM_CLOCK(start_count, count_rate)

      betanew = ( 2.0d0 / nc ) * beta

!     the ihat do loop, is to loop over all the possible Euclidean directions.
!     The imask do loop, is to ensure that the entire lattice is considered.

      do ihat=1,mu
         do imask=1,umask

!     calculate the staple in the ihat direction

            call staples(ur,ui,stapler,staplei,ihat,.true.,itype,uzero)

!
!  Insert loop here to loop over the three su2 subgroups twice
!
            do iloop=1,nloop

!     The Wilson action can be written as S(U)=Sum_p(Re(U*U_p))+constant
!                                             =Re(Tr(U*R))+constant
!     R is just the sum over the six plaquette namely the staples: stapler,staplei.

!     now we need to calculate the product of the staples time the link variable in
!     each direction U*R=ltsr+iltsi=( ur + iui ) * ( stapler + istaplei )
!                                  =(ur*stapler - ui*staplei)+i(ur*staplei + ui*staplei)
            do ic=1,nc
               do jc=1,nc
!
                  ltsr(:,:,:,:,ic,jc) =
     &                ur(:,:,:,:,ihat,ic,1) * stapler(:,:,:,:,1,jc) -
     &                ui(:,:,:,:,ihat,ic,1) * staplei(:,:,:,:,1,jc) +
     &                ur(:,:,:,:,ihat,ic,2) * stapler(:,:,:,:,2,jc) -
     &                ui(:,:,:,:,ihat,ic,2) * staplei(:,:,:,:,2,jc) +
     &                ur(:,:,:,:,ihat,ic,3) * stapler(:,:,:,:,3,jc) -
     &                ui(:,:,:,:,ihat,ic,3) * staplei(:,:,:,:,3,jc)
                  ltsi(:,:,:,:,ic,jc) =
     &                ur(:,:,:,:,ihat,ic,1) * staplei(:,:,:,:,1,jc) +
     &                ui(:,:,:,:,ihat,ic,1) * stapler(:,:,:,:,1,jc) +
     &                ur(:,:,:,:,ihat,ic,2) * staplei(:,:,:,:,2,jc) +
     &                ui(:,:,:,:,ihat,ic,2) * stapler(:,:,:,:,2,jc) +
     &                ur(:,:,:,:,ihat,ic,3) * staplei(:,:,:,:,3,jc) +
     &                ui(:,:,:,:,ihat,ic,3) * stapler(:,:,:,:,3,jc)
!
               end do
            end do

            phbsr(:,:,:,:,1,1) =
     &           ( ltsr(:,:,:,:,1,1) + ltsr(:,:,:,:,2,2) ) / 2.0d0
            phbsr(:,:,:,:,2,2) =   phbsr(:,:,:,:,1,1)
            phbsr(:,:,:,:,1,2) =
     &           ( ltsr(:,:,:,:,1,2) - ltsr(:,:,:,:,2,1) ) / 2.0d0
            phbsr(:,:,:,:,2,1) = - phbsr(:,:,:,:,1,2)

            phbsi(:,:,:,:,1,1) =
     &           ( ltsi(:,:,:,:,1,1) - ltsi(:,:,:,:,2,2) ) / 2.0d0
            phbsi(:,:,:,:,1,2) =
     &           ( ltsi(:,:,:,:,1,2) + ltsi(:,:,:,:,2,1) ) / 2.0d0
            phbsi(:,:,:,:,2,1) =   phbsi(:,:,:,:,1,2)
            phbsi(:,:,:,:,2,2) = - phbsi(:,:,:,:,1,1)

!     now calling the pseudo-heatbath algorithm to update, where the
!     mask is .true., the full QCD configuration

!     we call pseudoheat to form 2 SU(3) matrices.
!     These two matrices have the form matrixa = [ [SU(2)] 0 ],  matrixb[ 1    0    ].
!                                                [    0    1 ],         [ 0 [SU(2)] ]
!     Where [SU(2)] is hotwired via ursu2,uisu2 and ursu2,uisu2.

            call pseudoheat(ursu2,uisu2,phbsr,phbsi,mask,imask,betanew,ihat)

!     We next make a product of the matrices in such a way that
!     link Uprmsu3 = urprmsu3+iuiprmsu3 = ( ar + iai ) * ( ur + iui )
!                              = ( ar*ur - ai*ui ) + i( ar*ui + ai*ur)
!
!     The same multiplies are required for the ltsr and ltsi
!
            do ic=1,nc-1
               do jc=1,nc
                  where( mask(:,:,:,:,ihat,imask) )
              urprmsu3(:,:,:,:,ic,jc) = ( ursu2(:,:,:,:,ic,1) * ur(:,:,:,:,ihat,1,jc) +
     &                                    ursu2(:,:,:,:,ic,2) * ur(:,:,:,:,ihat,2,jc) -
     &                                    uisu2(:,:,:,:,ic,1) * ui(:,:,:,:,ihat,1,jc) -
     &                                    uisu2(:,:,:,:,ic,2) * ui(:,:,:,:,ihat,2,jc) )
              uiprmsu3(:,:,:,:,ic,jc) = ( ursu2(:,:,:,:,ic,1) * ui(:,:,:,:,ihat,1,jc) +
     &                                    ursu2(:,:,:,:,ic,2) * ui(:,:,:,:,ihat,2,jc) +
     &                                    uisu2(:,:,:,:,ic,1) * ur(:,:,:,:,ihat,1,jc) +
     &                                    uisu2(:,:,:,:,ic,2) * ur(:,:,:,:,ihat,2,jc) )
               ltsrprm(:,:,:,:,ic,jc) = ( ursu2(:,:,:,:,ic,1) * ltsr(:,:,:,:,1,jc) +
     &                                    ursu2(:,:,:,:,ic,2) * ltsr(:,:,:,:,2,jc) -
     &                                    uisu2(:,:,:,:,ic,1) * ltsi(:,:,:,:,1,jc) -
     &                                    uisu2(:,:,:,:,ic,2) * ltsi(:,:,:,:,2,jc) )
               ltsiprm(:,:,:,:,ic,jc) = ( ursu2(:,:,:,:,ic,1) * ltsi(:,:,:,:,1,jc) +
     &                                    ursu2(:,:,:,:,ic,2) * ltsi(:,:,:,:,2,jc) +
     &                                    uisu2(:,:,:,:,ic,1) * ltsr(:,:,:,:,1,jc) +
     &                                    uisu2(:,:,:,:,ic,2) * ltsr(:,:,:,:,2,jc) )
                  elsewhere
                     urprmsu3(:,:,:,:,ic,jc) = ur(:,:,:,:,ihat,ic,jc)
                     uiprmsu3(:,:,:,:,ic,jc) = ui(:,:,:,:,ihat,ic,jc)
                      ltsrprm(:,:,:,:,ic,jc) = ltsr(:,:,:,:,ic,jc)
                      ltsiprm(:,:,:,:,ic,jc) = ltsi(:,:,:,:,ic,jc)
                  end where
               end do
            end do

            do jc=1,nc
               urprmsu3(:,:,:,:,3,jc) = ur(:,:,:,:,ihat,3,jc)
               uiprmsu3(:,:,:,:,3,jc) = ui(:,:,:,:,ihat,3,jc)
                ltsrprm(:,:,:,:,3,jc) = ltsr(:,:,:,:,3,jc)
                ltsiprm(:,:,:,:,3,jc) = ltsi(:,:,:,:,3,jc)
            end do


            phbsr(:,:,:,:,1,1) =
     &           ( ltsrprm(:,:,:,:,2,2) + ltsrprm(:,:,:,:,3,3) ) / 2.0d0
            phbsr(:,:,:,:,2,2) =   phbsr(:,:,:,:,1,1)
            phbsr(:,:,:,:,1,2) =
     &           ( ltsrprm(:,:,:,:,2,3) - ltsrprm(:,:,:,:,3,2) ) / 2.0d0
            phbsr(:,:,:,:,2,1) = - phbsr(:,:,:,:,1,2)


            phbsi(:,:,:,:,1,1) =
     &           ( ltsiprm(:,:,:,:,2,2) - ltsiprm(:,:,:,:,3,3) ) / 2.0d0
            phbsi(:,:,:,:,1,2) =
     &           ( ltsiprm(:,:,:,:,2,3) + ltsiprm(:,:,:,:,3,2) ) / 2.0d0
            phbsi(:,:,:,:,2,1) =   phbsi(:,:,:,:,1,2)
            phbsi(:,:,:,:,2,2) = - phbsi(:,:,:,:,1,1)

!     Calling pseudoheat to obtain two new SU(2) matrices with the linkprmsu3*staples
!     calculated above.

            call pseudoheat(ursu2,uisu2,phbsr,phbsi,mask,imask,betanew,ihat)

!     Now mapping these variables to the full QCD link, ur and ui. For each
!     direction mu, the main loop imu.

            do ic=2,nc
               do jc=1,nc
                  where( mask(:,:,:,:,ihat,imask) )
              ur(:,:,:,:,ihat,ic,jc) = ( ursu2(:,:,:,:,ic-1,1) * urprmsu3(:,:,:,:,2,jc) +
     &                                   ursu2(:,:,:,:,ic-1,2) * urprmsu3(:,:,:,:,3,jc) -
     &                                   uisu2(:,:,:,:,ic-1,1) * uiprmsu3(:,:,:,:,2,jc) -
     &                                   uisu2(:,:,:,:,ic-1,2) * uiprmsu3(:,:,:,:,3,jc) )
              ui(:,:,:,:,ihat,ic,jc) = ( ursu2(:,:,:,:,ic-1,1) * uiprmsu3(:,:,:,:,2,jc) +
     &                                   ursu2(:,:,:,:,ic-1,2) * uiprmsu3(:,:,:,:,3,jc) +
     &                                   uisu2(:,:,:,:,ic-1,1) * urprmsu3(:,:,:,:,2,jc) +
     &                                   uisu2(:,:,:,:,ic-1,2) * urprmsu3(:,:,:,:,3,jc) )
                 ltsr(:,:,:,:,ic,jc) = ( ursu2(:,:,:,:,ic-1,1) * ltsrprm(:,:,:,:,2,jc) +
     &                                   ursu2(:,:,:,:,ic-1,2) * ltsrprm(:,:,:,:,3,jc) -
     &                                   uisu2(:,:,:,:,ic-1,1) * ltsiprm(:,:,:,:,2,jc) -
     &                                   uisu2(:,:,:,:,ic-1,2) * ltsiprm(:,:,:,:,3,jc) )
                 ltsi(:,:,:,:,ic,jc) = ( ursu2(:,:,:,:,ic-1,1) * ltsiprm(:,:,:,:,2,jc) +
     &                                   ursu2(:,:,:,:,ic-1,2) * ltsiprm(:,:,:,:,3,jc) +
     &                                   uisu2(:,:,:,:,ic-1,1) * ltsrprm(:,:,:,:,2,jc) +
     &                                   uisu2(:,:,:,:,ic-1,2) * ltsrprm(:,:,:,:,3,jc) )
                  elsewhere
                     ur(:,:,:,:,ihat,ic,jc) = urprmsu3(:,:,:,:,ic,jc)
                     ui(:,:,:,:,ihat,ic,jc) = uiprmsu3(:,:,:,:,ic,jc)
                        ltsr(:,:,:,:,ic,jc) = ltsrprm(:,:,:,:,ic,jc)
                        ltsi(:,:,:,:,ic,jc) = ltsiprm(:,:,:,:,ic,jc)
                  end where
               end do
            end do

            do jc=1,nc
               ur(:,:,:,:,ihat,1,jc) = urprmsu3(:,:,:,:,1,jc)
               ui(:,:,:,:,ihat,1,jc) = uiprmsu3(:,:,:,:,1,jc)
                  ltsr(:,:,:,:,1,jc) = ltsrprm(:,:,:,:,1,jc)
                  ltsi(:,:,:,:,1,jc) = ltsiprm(:,:,:,:,1,jc)
            end do

            phbsr(:,:,:,:,1,1) =
     &           ( ltsr(:,:,:,:,1,1) + ltsr(:,:,:,:,3,3) ) / 2.0d0
            phbsr(:,:,:,:,2,2) =   phbsr(:,:,:,:,1,1)
            phbsr(:,:,:,:,1,2) =
     &           ( ltsr(:,:,:,:,1,3) - ltsr(:,:,:,:,3,1) ) / 2.0d0
            phbsr(:,:,:,:,2,1) = - phbsr(:,:,:,:,1,2)

            phbsi(:,:,:,:,1,1) =
     &           ( ltsi(:,:,:,:,1,1) - ltsi(:,:,:,:,3,3) ) / 2.0d0
            phbsi(:,:,:,:,1,2) =
     &           ( ltsi(:,:,:,:,1,3) + ltsi(:,:,:,:,3,1) ) / 2.0d0
            phbsi(:,:,:,:,2,1) =   phbsi(:,:,:,:,1,2)
            phbsi(:,:,:,:,2,2) = - phbsi(:,:,:,:,1,1)


            call pseudoheat(ursu2,uisu2,phbsr,phbsi,mask,imask,betanew,ihat)

!     We next make a product of the matrices in such a way that
!     link Uprmsu3 = urprmsu3+iuiprmsu3 = ( ar + iai ) * ( ur + iui )
!                              = ( ar*ur - ai*ui ) + i( ar*ui + ai*ur)

            do ic=1,nc-1
               do jc=1,nc
                  ic3 = ic
                  if(ic3 == 2 ) ic3 = 3
                  where( mask(:,:,:,:,ihat,imask) )
            urprmsu3(:,:,:,:,ic3,jc) = ( ursu2(:,:,:,:,ic,1) * ur(:,:,:,:,ihat,1,jc) +
     &                                   ursu2(:,:,:,:,ic,2) * ur(:,:,:,:,ihat,3,jc) -
     &                                   uisu2(:,:,:,:,ic,1) * ui(:,:,:,:,ihat,1,jc) -
     &                                   uisu2(:,:,:,:,ic,2) * ui(:,:,:,:,ihat,3,jc) )
            uiprmsu3(:,:,:,:,ic3,jc) = ( ursu2(:,:,:,:,ic,1) * ui(:,:,:,:,ihat,1,jc) +
     &                                   ursu2(:,:,:,:,ic,2) * ui(:,:,:,:,ihat,3,jc) +
     &                                   uisu2(:,:,:,:,ic,1) * ur(:,:,:,:,ihat,1,jc) +
     &                                   uisu2(:,:,:,:,ic,2) * ur(:,:,:,:,ihat,3,jc) )
             ltsrprm(:,:,:,:,ic3,jc) = ( ursu2(:,:,:,:,ic,1) * ltsr(:,:,:,:,1,jc) +
     &                                   ursu2(:,:,:,:,ic,2) * ltsr(:,:,:,:,3,jc) -
     &                                   uisu2(:,:,:,:,ic,1) * ltsi(:,:,:,:,1,jc) -
     &                                   uisu2(:,:,:,:,ic,2) * ltsi(:,:,:,:,3,jc) )
             ltsiprm(:,:,:,:,ic3,jc) = ( ursu2(:,:,:,:,ic,1) * ltsi(:,:,:,:,1,jc) +
     &                                   ursu2(:,:,:,:,ic,2) * ltsi(:,:,:,:,3,jc) +
     &                                   uisu2(:,:,:,:,ic,1) * ltsr(:,:,:,:,1,jc) +
     &                                   uisu2(:,:,:,:,ic,2) * ltsr(:,:,:,:,3,jc) )
                  elsewhere
                     urprmsu3(:,:,:,:,ic3,jc) = ur(:,:,:,:,ihat,ic3,jc)
                     uiprmsu3(:,:,:,:,ic3,jc) = ui(:,:,:,:,ihat,ic3,jc)
                      ltsrprm(:,:,:,:,ic3,jc) = ltsr(:,:,:,:,ic3,jc)
                      ltsiprm(:,:,:,:,ic3,jc) = ltsi(:,:,:,:,ic3,jc)
                  end where
               end do
            end do

            do jc=1,nc
               ur(:,:,:,:,ihat,1,jc) = urprmsu3(:,:,:,:,1,jc)
               ui(:,:,:,:,ihat,1,jc) = uiprmsu3(:,:,:,:,1,jc)
               ur(:,:,:,:,ihat,3,jc) = urprmsu3(:,:,:,:,3,jc)
               ui(:,:,:,:,ihat,3,jc) = uiprmsu3(:,:,:,:,3,jc)
                  ltsr(:,:,:,:,1,jc) = ltsrprm(:,:,:,:,1,jc)
                  ltsi(:,:,:,:,1,jc) = ltsiprm(:,:,:,:,1,jc)
                  ltsr(:,:,:,:,3,jc) = ltsrprm(:,:,:,:,3,jc)
                  ltsi(:,:,:,:,3,jc) = ltsiprm(:,:,:,:,3,jc)
            end do

            end do
         end do
      end do

      CALL SYSTEM_CLOCK(end_count)
      elapsed_time = REAL(end_count - start_count) / REAL(count_rate)
!
!  Excluding assignments
!
      if (elapsed_time .ne. 0.0) then
        write(*,'(/,a,f7.3,a,/)') 'The elapsed time is',elapsed_time,
     &    ' seconds for one call to pseudosweep.'
      end if

      return

      end subroutine pseudosweep

      END MODULE GS_PSEUDOSWEEP
