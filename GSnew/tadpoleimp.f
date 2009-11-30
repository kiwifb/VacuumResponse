!
!     ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!     tadpoleimp subroutine will calculate the u_0 which is just the fourth
!     root of the average square plaquette. This uzero is uzero for the tadpole
!     improvement.
!
!     Authors: Frederi! D.R. Bonnet & Derek B. Leinweber
!     Date:    August 1998.
!
      MODULE GS_TADPOLEIMP

      CONTAINS

      subroutine tadpoleimp(ur,ui,uzero)

      USE GS_LATTICESIZE
      USE GS_SQUARES

      implicit none

!     global variables

      integer,parameter                                       :: nc=3               !sigma,color
      integer,parameter                                       :: mu=4               !direction

      double precision                                        :: uzero

      double precision,dimension(nx,ny,nz,nt,mu,nc,nc)        :: ur,ui
!HPF$ DISTRIBUTE ur(*,*,BLOCK,BLOCK,*,*,*)
!HPF$ DISTRIBUTE ui(*,*,BLOCK,BLOCK,*,*,*)

!     local variables

      double precision,dimension(nx,ny,nz,nt)                 :: action
!HPF$ DISTRIBUTE action(*,*,BLOCK,BLOCK)
      double precision,dimension(nx,ny,nz,nt,nc,nc)           :: squarer,squarei
!HPF$ DISTRIBUTE squarer(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE squarei(*,*,BLOCK,BLOCK,*,*)

      integer                                                 :: ihat,ic,kc

!     In the following do loop we are calculating the average plaquette.
!     when calling staples only three plaquette are calculated.
!     The ones in the positive direction ( not in the negative direction
!     because of the local=.false. logical condition).
!     at an arbitrary point we have
!     at ihat=1, plaquette 12,13,14. At ihat=2, plaquette 21,23,24
!     at ihat=3, plaquette 31,32,34. At ihat=4, plaquette 41,42,43

      action = 0.0d0
      uzero = 0.0d0

      do ihat=1,mu
         call squares(ur,ui,squarer,squarei,ihat,.false.)
         do ic=1,nc
            do kc=1,nc
               action(:,:,:,:) = action(:,:,:,:) +
     &              ( ur(:,:,:,:,ihat,ic,kc) * squarer(:,:,:,:,kc,ic) -
     &                ui(:,:,:,:,ihat,ic,kc) * squarei(:,:,:,:,kc,ic) )
            end do
         end do
      end do

!     now calculating the the fourth root of the plaquette. This uzero
!     makes up the tadpole improvement 1/3ReTr(square) consists of only
!     tadpole contributions. A large contribution of the tadpole contribution
!     can be eliminated by divinding every link operator by u_0

      uzero = ( sum(action) / ( nx*ny*nz*nt*nc*mu*(mu-1) ) ) ** ( 0.25d0 )

      return

      end subroutine tadpoleimp

      END MODULE GS_TADPOLEIMP
