!     subroutine that fixes the su3 links, using a smearing and cooling
!     techniques.
!
      MODULE GS_NEWFIXSU3

      CONTAINS

      subroutine Newfixsu3(ur,ui,urprm,uiprm,nsub)

      USE GS_LATTICESIZE
      USE GS_COOLING

      implicit none

!     global variables.

      integer,parameter                                       :: nc=3,ncsu2=2
      integer,parameter                                       :: mu=4

      integer                                                 :: nsub               !total # su2sub=nsub*3

      double precision,dimension(nx,ny,nz,nt,mu,nc,nc)        ::  urprm,uiprm
!HPF$ DISTRIBUTE urprm(*,*,BLOCK,BLOCK,*,*,*)
!HPF$ DISTRIBUTE uiprm(*,*,BLOCK,BLOCK,*,*,*)
      double precision,dimension(nx,ny,nz,nt,mu,nc,nc)        :: ur,ui
!HPF$ DISTRIBUTE ur(*,*,BLOCK,BLOCK,*,*,*)
!HPF$ DISTRIBUTE ui(*,*,BLOCK,BLOCK,*,*,*)

!     local variables.

      double precision,dimension(nx,ny,nz,nt,nc,nc)           ::stapler,staplei     !smeared links
!HPF$ DISTRIBUTE stapler(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE staplei(*,*,BLOCK,BLOCK,*,*)

      double precision,dimension(nx,ny,nz,nt,nc,nc)           :: ltsr,ltsi
      double precision,dimension(nx,ny,nz,nt,nc,nc)           :: urprmsu3,uiprmsu3
!HPF$ DISTRIBUTE ltsr(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE ltsi(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE urprmsu3(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE uiprmsu3(*,*,BLOCK,BLOCK,*,*)
      double precision,dimension(nx,ny,nz,nt,ncsu2,ncsu2)     :: phbsr,phbsi
      double precision,dimension(nx,ny,nz,nt,ncsu2,ncsu2)     :: ursu2,uisu2
!HPF$ DISTRIBUTE phbsr(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE phbsi(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE ursu2(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE uisu2(*,*,BLOCK,BLOCK,*,*)

      integer                     	                      :: ihat,isub
      integer                                                 :: of1
      integer                 		                      :: ic,jc,kc,ic3       !counters

!     starting the execution commands

!     first calculate the uprimes.

      do ihat=1,mu

!     now setting the urprm and uiprm dagger to the staples.

         do ic=1,nc
            do jc=1,nc
               stapler(:,:,:,:,ic,jc) =   urprm(:,:,:,:,ihat,jc,ic)
               staplei(:,:,:,:,ic,jc) = - uiprm(:,:,:,:,ihat,jc,ic)
            end do
         end do

!     now using the pseudosweep routine to cool the imformation at the
!     su(2) level.

!     here looping over the su(2) subgroups.

         do isub=1,nsub
!
!     case 1. a_1
!
!     The Wilson action can be written as S(U)=Sum_p(Re(U*U_p))+constant
!                                             =Re(Tr(U*R))+constant
!     R is just the sum over the six plaquette namely the staples: stapler,staplei.

            ltsr = 0.0d0
            ltsi = 0.0d0
            do ic=1,nc
               do jc=1,nc
                  do kc=1,nc
                     ltsr(:,:,:,:,ic,jc) = ltsr(:,:,:,:,ic,jc) +
     &                   ( ur(:,:,:,:,ihat,ic,kc) * stapler(:,:,:,:,kc,jc) -
     &                     ui(:,:,:,:,ihat,ic,kc) * staplei(:,:,:,:,kc,jc) )
                     ltsi(:,:,:,:,ic,jc) = ltsi(:,:,:,:,ic,jc) +
     &                   ( ur(:,:,:,:,ihat,ic,kc) * staplei(:,:,:,:,kc,jc) +
     &                     ui(:,:,:,:,ihat,ic,kc) * stapler(:,:,:,:,kc,jc) )
                  end do
               end do
            end do

            of1 = 0

            phbsr(:,:,:,:,1,1) =
     &           ( ltsr(:,:,:,:,1+of1,1+of1) + ltsr(:,:,:,:,2+of1,2+of1) ) / 2.0d0
            phbsr(:,:,:,:,2,2) =   phbsr(:,:,:,:,1,1)
            phbsr(:,:,:,:,1,2) =
     &           ( ltsr(:,:,:,:,1+of1,2+of1) - ltsr(:,:,:,:,2+of1,1+of1) ) / 2.0d0
            phbsr(:,:,:,:,2,1) = - phbsr(:,:,:,:,1,2)

            phbsi(:,:,:,:,1,1) =
     &           ( ltsi(:,:,:,:,1+of1,1+of1) - ltsi(:,:,:,:,2+of1,2+of1) ) / 2.0d0
            phbsi(:,:,:,:,1,2) =
     &           ( ltsi(:,:,:,:,1+of1,2+of1) + ltsi(:,:,:,:,2+of1,1+of1) ) / 2.0d0
            phbsi(:,:,:,:,2,1) =   phbsi(:,:,:,:,1,2)
            phbsi(:,:,:,:,2,2) = - phbsi(:,:,:,:,1,1)

            call cooling(ursu2,uisu2,phbsr,phbsi)

            do ic=1,nc-1
               do jc=1,nc
                  urprmsu3(:,:,:,:,ic,jc) = ( ursu2(:,:,:,:,ic,1) * ur(:,:,:,:,ihat,1,jc) +
     &                                        ursu2(:,:,:,:,ic,2) * ur(:,:,:,:,ihat,2,jc) -
     &                                        uisu2(:,:,:,:,ic,1) * ui(:,:,:,:,ihat,1,jc) -
     &                                        uisu2(:,:,:,:,ic,2) * ui(:,:,:,:,ihat,2,jc) )
                  uiprmsu3(:,:,:,:,ic,jc) = ( ursu2(:,:,:,:,ic,1) * ui(:,:,:,:,ihat,1,jc) +
     &                                        ursu2(:,:,:,:,ic,2) * ui(:,:,:,:,ihat,2,jc) +
     &                                        uisu2(:,:,:,:,ic,1) * ur(:,:,:,:,ihat,1,jc) +
     &                                        uisu2(:,:,:,:,ic,2) * ur(:,:,:,:,ihat,2,jc) )
               end do
            end do

            do jc=1,nc
               urprmsu3(:,:,:,:,3,jc) = ur(:,:,:,:,ihat,3,jc)
               uiprmsu3(:,:,:,:,3,jc) = ui(:,:,:,:,ihat,3,jc)
            end do
!
!     case 2. a_2
!
!     here we calculate the next bit of the product, link*staples.
!     We use the old staples to multiply it with Uprmsu3 = a_1 * U
!     lts = matrixa = [ [SU(2)] 0 ]*lts[ ele in SU(3)    0    ].
!                     [    0    1 ]    [      0          0    ]

            ltsr = 0.0d0
            ltsi = 0.0d0
            do ic=1,nc
               do jc=1,nc
                  do kc=1,nc
                     ltsr(:,:,:,:,ic,jc) = ltsr(:,:,:,:,ic,jc) +
     &                   ( urprmsu3(:,:,:,:,ic,kc) * stapler(:,:,:,:,kc,jc) -
     &                     uiprmsu3(:,:,:,:,ic,kc) * staplei(:,:,:,:,kc,jc) )
                     ltsi(:,:,:,:,ic,jc) = ltsi(:,:,:,:,ic,jc) +
     &                   ( urprmsu3(:,:,:,:,ic,kc) * staplei(:,:,:,:,kc,jc) +
     &                     uiprmsu3(:,:,:,:,ic,kc) * stapler(:,:,:,:,kc,jc) )
                  end do
               end do
            end do

            of1 = 1

            phbsr(:,:,:,:,1,1) =
     &           ( ltsr(:,:,:,:,1+of1,1+of1) + ltsr(:,:,:,:,2+of1,2+of1) ) / 2.0d0
            phbsr(:,:,:,:,2,2) =   phbsr(:,:,:,:,1,1)
            phbsr(:,:,:,:,1,2) =
     &           ( ltsr(:,:,:,:,1+of1,2+of1) - ltsr(:,:,:,:,2+of1,1+of1) ) / 2.0d0
            phbsr(:,:,:,:,2,1) = - phbsr(:,:,:,:,1,2)


            phbsi(:,:,:,:,1,1) =
     &           ( ltsi(:,:,:,:,1+of1,1+of1) - ltsi(:,:,:,:,2+of1,2+of1) ) / 2.0d0
            phbsi(:,:,:,:,1,2) =
     &           ( ltsi(:,:,:,:,1+of1,2+of1) + ltsi(:,:,:,:,2+of1,1+of1) ) / 2.0d0
            phbsi(:,:,:,:,2,1) =   phbsi(:,:,:,:,1,2)
            phbsi(:,:,:,:,2,2) = - phbsi(:,:,:,:,1,1)

            call cooling(ursu2,uisu2,phbsr,phbsi)

            do ic=2,nc
               do jc=1,nc
                  ur(:,:,:,:,ihat,ic,jc) = ( ursu2(:,:,:,:,ic-1,1) * urprmsu3(:,:,:,:,2,jc) +
     &                                       ursu2(:,:,:,:,ic-1,2) * urprmsu3(:,:,:,:,3,jc) -
     &                                       uisu2(:,:,:,:,ic-1,1) * uiprmsu3(:,:,:,:,2,jc) -
     &                                       uisu2(:,:,:,:,ic-1,2) * uiprmsu3(:,:,:,:,3,jc) )
                  ui(:,:,:,:,ihat,ic,jc) = ( ursu2(:,:,:,:,ic-1,1) * uiprmsu3(:,:,:,:,2,jc) +
     &                                       ursu2(:,:,:,:,ic-1,2) * uiprmsu3(:,:,:,:,3,jc) +
     &                                       uisu2(:,:,:,:,ic-1,1) * urprmsu3(:,:,:,:,2,jc) +
     &                                       uisu2(:,:,:,:,ic-1,2) * urprmsu3(:,:,:,:,3,jc) )
               end do
            end do

            do jc=1,nc
               ur(:,:,:,:,ihat,1,jc) = urprmsu3(:,:,:,:,1,jc)
               ui(:,:,:,:,ihat,1,jc) = uiprmsu3(:,:,:,:,1,jc)
            end do
!
!     case 3. a_3
!
!     forming another SU(2) subgroup of the form
!     These two matrices have the form matrixa = [ x  0   x ]
!                                                [ 0  1   0 ]
!                                                [ x  0   x ]
            ltsr = 0.0d0
            ltsi = 0.0d0
            do ic=1,nc
               do jc=1,nc
                  do kc=1,nc
                     ltsr(:,:,:,:,ic,jc) = ltsr(:,:,:,:,ic,jc) +
     &                   ( ur(:,:,:,:,ihat,ic,kc) * stapler(:,:,:,:,kc,jc) -
     &                     ui(:,:,:,:,ihat,ic,kc) * staplei(:,:,:,:,kc,jc) )
                     ltsi(:,:,:,:,ic,jc) = ltsi(:,:,:,:,ic,jc) +
     &                   ( ur(:,:,:,:,ihat,ic,kc) * staplei(:,:,:,:,kc,jc) +
     &                     ui(:,:,:,:,ihat,ic,kc) * stapler(:,:,:,:,kc,jc) )
                  end do
               end do
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

            call cooling(ursu2,uisu2,phbsr,phbsi)

!     We next make a product of the matrices in such a way that
!     link Uprmsu3 = urprmsu3+iuiprmsu3 = ( ar + iai ) * ( ur + iui )
!                              = ( ar*ur - ai*ui ) + i( ar*ui + ai*ur)
!     by hotwiring the matrix indices for optimization.

            do ic=1,nc-1
               do jc=1,nc
                  ic3 = ic
                  if(ic3 == 2 ) ic3 = 3
                  urprmsu3(:,:,:,:,ic3,jc) = ( ursu2(:,:,:,:,ic,1) * ur(:,:,:,:,ihat,1,jc) +
     &                                         ursu2(:,:,:,:,ic,2) * ur(:,:,:,:,ihat,3,jc) -
     &                                         uisu2(:,:,:,:,ic,1) * ui(:,:,:,:,ihat,1,jc) -
     &                                         uisu2(:,:,:,:,ic,2) * ui(:,:,:,:,ihat,3,jc) )
                  uiprmsu3(:,:,:,:,ic3,jc) = ( ursu2(:,:,:,:,ic,1) * ui(:,:,:,:,ihat,1,jc) +
     &                                         ursu2(:,:,:,:,ic,2) * ui(:,:,:,:,ihat,3,jc) +
     &                                         uisu2(:,:,:,:,ic,1) * ur(:,:,:,:,ihat,1,jc) +
     &                                         uisu2(:,:,:,:,ic,2) * ur(:,:,:,:,ihat,3,jc) )
               end do
            end do

            do jc=1,nc
               ur(:,:,:,:,ihat,1,jc) = urprmsu3(:,:,:,:,1,jc)
               ui(:,:,:,:,ihat,1,jc) = uiprmsu3(:,:,:,:,1,jc)
               ur(:,:,:,:,ihat,3,jc) = urprmsu3(:,:,:,:,3,jc)
               ui(:,:,:,:,ihat,3,jc) = uiprmsu3(:,:,:,:,3,jc)
            end do

         end do

      end do

      return

      end subroutine Newfixsu3

      END MODULE GS_NEWFIXSU3
