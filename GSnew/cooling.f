!
!     cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!     subroutine that implements the pseudo-heatbath algorithm
!     Author: Frederi! D.R. Bonnet & D.B. Leinweber
!     Date: February 1999
!
!
      MODULE GS_COOLING

      CONTAINS

      subroutine cooling(urnewsu2,uinewsu2,phbsr,phbsi)

      USE GS_LATTICESIZE

      implicit none

!     global variables

      integer,parameter                                       :: ncsu2=2

      double precision,dimension(nx,ny,nz,nt,ncsu2,ncsu2)     :: urnewsu2,uinewsu2
      double precision,dimension(nx,ny,nz,nt,ncsu2,ncsu2)     :: phbsr,phbsi
!HPF$ DISTRIBUTE urnewsu2(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE uinewsu2(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE phbsr(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE phbsi(*,*,BLOCK,BLOCK,*,*)

!     local variables

      double precision,dimension(nx,ny,nz,nt)                 :: k
!HPF$ DISTRIBUTE k(*,*,BLOCK,BLOCK)

      integer                                                 :: ic,jc

!     start of the execution commands

!     calculates the determinant k=|\sum_{\alpha=1}^6\widetilde{U}_{\alpha}|^{1/2}
!     where $\widetilde{U}_{\alpha}\equiv$ the six product of the three links
!     variable which interact with the link in question, i.e. stapler and staplei

        k = sqrt( abs( phbsr(:,:,:,:,1,1) * phbsr(:,:,:,:,2,2) -
     &                 phbsr(:,:,:,:,1,2) * phbsr(:,:,:,:,2,1) -
     &                 phbsi(:,:,:,:,1,1) * phbsi(:,:,:,:,2,2) +
     &                 phbsi(:,:,:,:,1,2) * phbsi(:,:,:,:,2,1) ) )

!     this calculates U --> U'= staple^{\dag} / k, the U (full link) coming out
!     of su2random is here replaced by U'(partial link, depends on ihat the direction)
!     this is for minimising the local action such that U*{\overline{U}}^{\dag}=I

      urnewsu2 = 0.0d0
      uinewsu2 = 0.0d0
      do ic=1,ncsu2
         do jc=1,ncsu2
                  urnewsu2(:,:,:,:,ic,jc) =   phbsr(:,:,:,:,jc,ic) / k(:,:,:,:)
                  uinewsu2(:,:,:,:,ic,jc) = - phbsi(:,:,:,:,jc,ic) / k(:,:,:,:)
         end do
      end do

      return

      end subroutine cooling

      END MODULE GS_COOLING
