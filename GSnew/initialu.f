!
!     cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!     initializes the link variables
!
!     Author: Frederic D.R. Bonnet and Derek Leinweber
!     Date: June 1998
!
      MODULE GS_INITIALU

      CONTAINS

      subroutine initialu(ur,ui)

      USE GS_LATTICESIZE

      implicit none

!     global variable

      integer,parameter                                       :: nc=3
      integer,parameter                                       :: mu=4

      double precision,dimension(nx,ny,nz,nt,mu,nc,nc)        :: ur,ui
!HPF$ DISTRIBUTE ur(*,*,BLOCK,BLOCK,*,*,*)
!HPF$ DISTRIBUTE ui(*,*,BLOCK,BLOCK,*,*,*)

!     local variables

      integer                                                 :: ic

      ur = 0.0d0
      ui = 0.0d0
      do ic = 1, nc
        ur(:,:,:,:,:,ic,ic) = 1.0d0
      end do

      return

      end subroutine initialu

      END MODULE GS_INITIALU
