!
!     cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!     subroutine su3random calculates the random su3 configuration
!     it gets a complete set of random number for the lattice in the SU(3)
!     group.
!
!     Authors: Frederi! D.R. Bonnet & Derek B. Leinweber
!     Date: July 1998.
!
      MODULE GS_SU3RANDOM

      CONTAINS

      subroutine su3random(ur,ui)

      USE GS_LATTICESIZE
      USE GS_SU2RANDOM

      implicit none

!     global variables

      integer,parameter                                       :: ncsu2=2,nc=3
      integer,parameter                                       :: mu=4

      double precision,dimension(nx,ny,nz,nt,mu,nc,nc)        :: ur,ui
!HPF$ DISTRIBUTE ur(*,*,BLOCK,BLOCK,*,*,*)
!HPF$ DISTRIBUTE ui(*,*,BLOCK,BLOCK,*,*,*)

!     local variables

      double precision,dimension(nx,ny,nz,nt,mu,nc,nc)        :: urprm,uiprm
!HPF$ DISTRIBUTE urprm(*,*,BLOCK,BLOCK,*,*,*)
!HPF$ DISTRIBUTE uiprm(*,*,BLOCK,BLOCK,*,*,*)
      double precision,dimension(nx,ny,nz,nt,mu,ncsu2,ncsu2)  :: ursu2,uisu2
!HPF$ DISTRIBUTE ursu2(*,*,BLOCK,BLOCK,*,*,*)
!HPF$ DISTRIBUTE uisu2(*,*,BLOCK,BLOCK,*,*,*)

      integer                                                 :: ic,jc

!     start of the executable commands

!     we call su2random to form a matrix called matrixa,matrixb.
!     these two matrix have the form matrixa = [ [SU(2)] 0 ],  matrixb[ 1    0    ].
!                                              [    0    1 ],         [ 0 [SU(2)] ]
!     Where [SU(2)] is hotwired vi ursu2,uisu2.
!     Both matrixa and matrixb have a double precision and imaginary part. matrixar,matrixai,
!     similarly for matrixb

      call su2random(ursu2,uisu2)

!     we next make a product of the matrices in such a way that
!     link Urpm=urprm+iuiprm=( ar + iai ) * ( ur + iui )
!                           =( ar*ur - ai*ui ) + i( ar*ui + ai*ur)
!     by hardwiring the matrix indices for optimization

      do ic=1,nc-1
        do jc=1,nc
          urprm(:,:,:,:,:,ic,jc) =
     &      ( ursu2(:,:,:,:,:,ic,1) * ur(:,:,:,:,:,1,jc) +
     &        ursu2(:,:,:,:,:,ic,2) * ur(:,:,:,:,:,2,jc) -
     &        uisu2(:,:,:,:,:,ic,1) * ui(:,:,:,:,:,1,jc) -
     &        uisu2(:,:,:,:,:,ic,2) * ui(:,:,:,:,:,2,jc) )
          uiprm(:,:,:,:,:,ic,jc) =
     &      ( ursu2(:,:,:,:,:,ic,1) * ui(:,:,:,:,:,1,jc) +
     &        ursu2(:,:,:,:,:,ic,2) * ui(:,:,:,:,:,2,jc) +
     &        uisu2(:,:,:,:,:,ic,1) * ur(:,:,:,:,:,1,jc) +
     &        uisu2(:,:,:,:,:,ic,2) * ur(:,:,:,:,:,2,jc) )
        end do
      end do

      do jc=1,nc
        urprm(:,:,:,:,:,3,jc) = ur(:,:,:,:,:,3,jc)
        uiprm(:,:,:,:,:,3,jc) = ui(:,:,:,:,:,3,jc)
      end do

!     we next make a product of the matrices in such a way that link
!     Udblerpm=urdbleprm+iuidbleprm = ( br + ibi ) * ( urprm + iuiprm )
!                                   = ( br*urprm - bi*uiprm ) + i( br*uiprm + bi*urprm)

      call su2random(ursu2,uisu2)

      do ic=2,nc
        do jc=1,nc
          ur(:,:,:,:,:,ic,jc) =
     &      ( ursu2(:,:,:,:,:,ic-1,1) * urprm(:,:,:,:,:,2,jc) +
     &        ursu2(:,:,:,:,:,ic-1,2) * urprm(:,:,:,:,:,3,jc) -
     &        uisu2(:,:,:,:,:,ic-1,1) * uiprm(:,:,:,:,:,2,jc) -
     &        uisu2(:,:,:,:,:,ic-1,2) * uiprm(:,:,:,:,:,3,jc) )
          ui(:,:,:,:,:,ic,jc) =
     &      ( ursu2(:,:,:,:,:,ic-1,1) * uiprm(:,:,:,:,:,2,jc) +
     &        ursu2(:,:,:,:,:,ic-1,2) * uiprm(:,:,:,:,:,3,jc) +
     &        uisu2(:,:,:,:,:,ic-1,1) * urprm(:,:,:,:,:,2,jc) +
     &        uisu2(:,:,:,:,:,ic-1,2) * urprm(:,:,:,:,:,3,jc) )
        end do
      end do

      do jc=1,nc
        ur(:,:,:,:,:,1,jc) = urprm(:,:,:,:,:,1,jc)
        ui(:,:,:,:,:,1,jc) = uiprm(:,:,:,:,:,1,jc)
      end do

      return

      end subroutine su3random

      END MODULE GS_SU3RANDOM
