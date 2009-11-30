!
!
!------------------------------------------------------------------------------------------
!
!     MaskWilson subroutine establishes the x y z and t masks.
!     This mask is only to be used for the use of ordinary plaquette action
!
!  Author: Frederi! Bonnet  fbonnet@physics.adelaide.edu.au
!          Updated 14 Sept. 2000  by DBL
!
      MODULE GS_MASKWILSON

      CONTAINS

      subroutine MaskWilson(mask)

      USE GS_LATTICESIZE

      implicit none
      integer,parameter                                       :: mu=4
      integer,parameter                                       :: nmask=16
      logical,dimension(nx,ny,nz,nt,mu,nmask)                 :: mask
!HPF$ DISTRIBUTE mask(*,*,BLOCK,BLOCK,*,*)

!     local variables

      integer                                                 :: ix,iy,iz,it
!
!------------------------------------------------------------------------------------------
!
      mask = .false.
!
!------------------------------------------------------------------------------------------
!
!  Our first idea
!
!      forall( ix=1:nx,iy=1:ny,iz=1:nz,it=1:nt,mod(iy+iz+it,2) .eq. 0 )
!     &     mask(ix,iy,iz,it,1,1)=.true.
!      forall( ix=1:nx,iy=1:ny,iz=1:nz,it=1:nt,mod(ix+iz+it,2) .eq. 0 )
!     &     mask(ix,iy,iz,it,2,1)=.true.
!      forall( ix=1:nx,iy=1:ny,iz=1:nz,it=1:nt,mod(ix+iy+it,2) .eq. 0 )
!     &     mask(ix,iy,iz,it,3,1)=.true.
!      forall( ix=1:nx,iy=1:ny,iz=1:nz,it=1:nt,mod(ix+iy+iz,2) .eq. 0 )
!     &     mask(ix,iy,iz,it,4,1)=.true.
!
!  Checker Board this time
!
      forall( ix=1:nx,iy=1:ny,iz=1:nz,it=1:nt,mod(ix+iy+iz+it,2) .eq. 0 )
     &     mask(ix,iy,iz,it,1,1)=.true.
      forall( ix=1:nx,iy=1:ny,iz=1:nz,it=1:nt,mod(ix+iy+iz+it,2) .eq. 0 )
     &     mask(ix,iy,iz,it,2,1)=.true.
      forall( ix=1:nx,iy=1:ny,iz=1:nz,it=1:nt,mod(ix+iy+iz+it,2) .eq. 0 )
     &     mask(ix,iy,iz,it,3,1)=.true.
      forall( ix=1:nx,iy=1:ny,iz=1:nz,it=1:nt,mod(ix+iy+iz+it,2) .eq. 0 )
     &     mask(ix,iy,iz,it,4,1)=.true.

      mask(:,:,:,:,:,2) = .not. mask(:,:,:,:,:,1)

      return

      end subroutine MaskWilson

      END MODULE GS_MASKWILSON
