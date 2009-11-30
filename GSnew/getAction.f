!
!     cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!     GetAction subroutine will calculate the total action, which is a sum
!     over all elementary action in the lattice.
!
!     Author: Frederi! D.R. Bonnet & D.B. Leinweber: date: March 1999.
!
      MODULE GS_GETACTION

      CONTAINS

      subroutine GetAction(ur,ui,action,plaqbar,itype,uzero)

      USE GS_LATTICESIZE
      USE GS_STAPLES

      implicit none

!     global variables

      integer,parameter                                       :: nc=3
      integer,parameter                                       :: mu=4

      integer                                                 :: itype
      double precision                                        :: uzero
      double precision                                        :: plaqbar

      double precision,dimension(nx,ny,nz,nt,mu,nc,nc)        :: ur,ui
!HPF$ DISTRIBUTE ur(*,*,BLOCK,BLOCK,*,*,*)
!HPF$ DISTRIBUTE ui(*,*,BLOCK,BLOCK,*,*,*)
      double precision,dimension(nx,ny,nz,nt)                 :: action
!HPF$ DISTRIBUTE action(*,*,BLOCK,BLOCK)

!     local variables

      double precision                                        :: factor

      double precision,dimension(nx,ny,nz,nt,nc,nc)           :: stapler,staplei
!HPF$ DISTRIBUTE stapler(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE staplei(*,*,BLOCK,BLOCK,*,*)

      integer                                                 :: ihat,ic,kc

!     In the following do loop we are calculating the average plaquette.
!     when calling staples only three plaquette are calculated.
!     The ones in the positive direction ( not in the negative direction
!     because of the local=.false. logical condition).
!     at an arbitrary point we have
!     at ihat=1, plaquette 12,13,14. At ihat=2, plaquette 21,23,24
!     at ihat=3, plaquette 31,32,34. At ihat=4, plaquette 41,42,43

!     Initialize action to zero here.

      if( itype == 0 ) factor = 1.0d0
      if( itype == 1 ) factor = 5.0d0 / 3.0d0 - 1.0d0 / (6.0d0 * uzero**2)

      action  = 0.0d0
      plaqbar = 0.0d0

      stapler = 0.0d0
      staplei = 0.0d0

      do ihat=1,mu
         call staples(ur,ui,stapler,staplei,ihat,.false.,itype,uzero)
         do ic=1,nc
            do kc=1,nc
               action(:,:,:,:)    = action(:,:,:,:)    +
     &              ( ur(:,:,:,:,ihat,ic,kc) * stapler(:,:,:,:,kc,ic) -
     &                ui(:,:,:,:,ihat,ic,kc) * staplei(:,:,:,:,kc,ic) )
            end do
         end do
      end do

!     each subroutine call to staples calculates 3 plaquettes (mu-1)
!     when called with .false., and there are 4 calls, one for each direction (mu).

      action = factor - action / ( nc * mu*(mu-1) )

      plaqbar = sum(action) / ( nx*ny*nz*nt )

      return

      end subroutine GetAction

      END MODULE GS_GETACTION
