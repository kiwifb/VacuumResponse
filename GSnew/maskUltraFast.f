!
!
!------------------------------------------------------------------------------------------
!
!     Subroutine MaskUltraFast modifies MaskUltra to nest the plaquette plus rectangle
!       action.  This is the mask of choice as only 4 masks per link direction are
!       required.
!
!     Authors: Frederi! Bonnet           fbonnet@physics.adelaide.edu.au
!              Sundance Bilson-Thompson  sbilson@physics.adelaide.edu.au
!              Derek B. Leinweber        dleinweb@physics.adelaide.edu.au
!              14 Sept. 2000
!
      MODULE GS_MASKULTRAFAST

      CONTAINS

      subroutine MaskUltraFast(mask)

      USE GS_LATTICESIZE

      implicit none
!     global variables

      integer,parameter                       :: mu=4
      integer,parameter                       :: nmask=16

      logical,dimension(nx,ny,nz,nt,mu,nmask) :: mask
!HPF$ DISTRIBUTE mask(*,*,BLOCK,BLOCK,*,*)

!     local variables

      integer,parameter                       :: ncube=4
      integer,parameter                       :: nlink=2
      integer,parameter                       :: nspace=3
      integer                                 :: ix,iy,iz,it
      integer                                 :: imask
      integer                                 :: idir
!
!------------------------------------------------------------------------------------------
!
      mask = .false.
!
!------------------------------------------------------------------------------------------
!
!

!     mask 1 for the x direction

      do ix=1,nx,nlink
         do iy=1,ny,ncube
            do iz=1,nz,ncube
               do it=1,nt,ncube
!     plane 1
                  mask( ix , iy+0 , iz+0 , it+0 , 1 , 1 ) = .true.
                  mask( ix , iy+1 , iz+1 , it+0 , 1 , 1 ) = .true.
                  mask( ix , iy+2 , iz+2 , it+0 , 1 , 1 ) = .true.
                  mask( ix , iy+3 , iz+3 , it+0 , 1 , 1 ) = .true.
!     plane 2
                  mask( ix , iy+0 , iz+1 , it+1 , 1 , 1 ) = .true.
                  mask( ix , iy+1 , iz+2 , it+1 , 1 , 1 ) = .true.
                  mask( ix , iy+2 , iz+3 , it+1 , 1 , 1 ) = .true.
                  mask( ix , iy+3 , iz+0 , it+1 , 1 , 1 ) = .true.
!     plane 3
                  mask( ix , iy+0 , iz+2 , it+2 , 1 , 1 ) = .true.
                  mask( ix , iy+1 , iz+3 , it+2 , 1 , 1 ) = .true.
                  mask( ix , iy+2 , iz+0 , it+2 , 1 , 1 ) = .true.
                  mask( ix , iy+3 , iz+1 , it+2 , 1 , 1 ) = .true.
!     plane 4
                  mask( ix , iy+0 , iz+3 , it+3 , 1 , 1 ) = .true.
                  mask( ix , iy+1 , iz+0 , it+3 , 1 , 1 ) = .true.
                  mask( ix , iy+2 , iz+1 , it+3 , 1 , 1 ) = .true.
                  mask( ix , iy+3 , iz+2 , it+3 , 1 , 1 ) = .true.

               end do
            end do
         end do
      end do

!
!  Now shift this mask 1 step in x direction and 2 in y, z, and t directions
!    and OR it to the first
!
      mask(:,:,:,:,1,1) = mask(:,:,:,:,1,1) .or.
     &         cshift(
     &           cshift(
     &             cshift(
     &               cshift(mask(:,:,:,:,1,1), dim=1, shift=1),
     &             dim=2, shift=2),
     &           dim=3, shift=2),
     &         dim=4, shift=2)


!     mask 1 for the y direction

      do ix=1,nx,ncube
         do iy=1,ny,nlink
            do iz=1,nz,ncube
               do it=1,nt,ncube
!     plane 1
                  mask( ix+0 , iy , iz+0 , it+0 , 2 , 1 ) = .true.
                  mask( ix+1 , iy , iz+1 , it+0 , 2 , 1 ) = .true.
                  mask( ix+2 , iy , iz+2 , it+0 , 2 , 1 ) = .true.
                  mask( ix+3 , iy , iz+3 , it+0 , 2 , 1 ) = .true.
!     plane 2
                  mask( ix+0 , iy , iz+1 , it+1 , 2 , 1 ) = .true.
                  mask( ix+1 , iy , iz+2 , it+1 , 2 , 1 ) = .true.
                  mask( ix+2 , iy , iz+3 , it+1 , 2 , 1 ) = .true.
                  mask( ix+3 , iy , iz+0 , it+1 , 2 , 1 ) = .true.
!     plane 3
                  mask( ix+0 , iy , iz+2 , it+2 , 2 , 1 ) = .true.
                  mask( ix+1 , iy , iz+3 , it+2 , 2 , 1 ) = .true.
                  mask( ix+2 , iy , iz+0 , it+2 , 2 , 1 ) = .true.
                  mask( ix+3 , iy , iz+1 , it+2 , 2 , 1 ) = .true.
!     plane 4
                  mask( ix+0 , iy , iz+3 , it+3 , 2 , 1 ) = .true.
                  mask( ix+1 , iy , iz+0 , it+3 , 2 , 1 ) = .true.
                  mask( ix+2 , iy , iz+1 , it+3 , 2 , 1 ) = .true.
                  mask( ix+3 , iy , iz+2 , it+3 , 2 , 1 ) = .true.

               end do
            end do
         end do
      end do

!
!  Now shift this mask 1 step in y direction and 2 in x, z, and t directions
!    and OR it to the first
!
      mask(:,:,:,:,2,1) = mask(:,:,:,:,2,1) .or.
     &         cshift(
     &           cshift(
     &             cshift(
     &               cshift(mask(:,:,:,:,2,1), dim=2, shift=1),
     &             dim=1, shift=2),
     &           dim=3, shift=2),
     &         dim=4, shift=2)


!     mask 1 for the z direction

      do ix=1,nx,ncube
         do iy=1,ny,ncube
            do iz=1,nz,nlink
               do it=1,nt,ncube
!     plane 1
                  mask( ix+0 , iy+0 , iz , it+0 , 3 , 1 ) = .true.
                  mask( ix+1 , iy+1 , iz , it+0 , 3 , 1 ) = .true.
                  mask( ix+2 , iy+2 , iz , it+0 , 3 , 1 ) = .true.
                  mask( ix+3 , iy+3 , iz , it+0 , 3 , 1 ) = .true.
!     plane 2
                  mask( ix+0 , iy+1 , iz , it+1 , 3 , 1 ) = .true.
                  mask( ix+1 , iy+2 , iz , it+1 , 3 , 1 ) = .true.
                  mask( ix+2 , iy+3 , iz , it+1 , 3 , 1 ) = .true.
                  mask( ix+3 , iy+0 , iz , it+1 , 3 , 1 ) = .true.
!     plane 3
                  mask( ix+0 , iy+2 , iz , it+2 , 3 , 1 ) = .true.
                  mask( ix+1 , iy+3 , iz , it+2 , 3 , 1 ) = .true.
                  mask( ix+2 , iy+0 , iz , it+2 , 3 , 1 ) = .true.
                  mask( ix+3 , iy+1 , iz , it+2 , 3 , 1 ) = .true.
!     plane 4
                  mask( ix+0 , iy+3 , iz , it+3 , 3 , 1 ) = .true.
                  mask( ix+1 , iy+0 , iz , it+3 , 3 , 1 ) = .true.
                  mask( ix+2 , iy+1 , iz , it+3 , 3 , 1 ) = .true.
                  mask( ix+3 , iy+2 , iz , it+3 , 3 , 1 ) = .true.
               end do
            end do
         end do
      end do

!
!  Now shift this mask 1 step in z direction and 2 in x, y, and t directions
!    and OR it to the first
!
      mask(:,:,:,:,3,1) = mask(:,:,:,:,3,1) .or.
     &         cshift(
     &           cshift(
     &             cshift(
     &               cshift(mask(:,:,:,:,3,1), dim=3, shift=1),
     &             dim=1, shift=2),
     &           dim=2, shift=2),
     &         dim=4, shift=2)

!     mask 2, 3, 4

      do idir = 1, nspace
        do imask = 2, ncube
          mask(:,:,:,:,idir,imask) = cshift( mask(:,:,:,:,idir,imask-1), dim=4, shift=1 )
        end do
      end do
!
!------------------------------------------------------------------------------------------
!
!     doing 3 steps in the t direction

!     mask 1

      do ix=1,nx,ncube
         do iy=1,ny,ncube
            do iz=1,nz,ncube
               do it=1,nt,nlink
!     plane 1
                  mask( ix+0 , iy+0 , iz+0 , it , 4 , 1 ) = .true.
                  mask( ix+0 , iy+1 , iz+1 , it , 4 , 1 ) = .true.
                  mask( ix+0 , iy+2 , iz+2 , it , 4 , 1 ) = .true.
                  mask( ix+0 , iy+3 , iz+3 , it , 4 , 1 ) = .true.
!     plane 2
                  mask( ix+1 , iy+0 , iz+1 , it , 4 , 1 ) = .true.
                  mask( ix+1 , iy+1 , iz+2 , it , 4 , 1 ) = .true.
                  mask( ix+1 , iy+2 , iz+3 , it , 4 , 1 ) = .true.
                  mask( ix+1 , iy+3 , iz+0 , it , 4 , 1 ) = .true.
!     plane 3
                  mask( ix+2 , iy+0 , iz+2 , it , 4 , 1 ) = .true.
                  mask( ix+2 , iy+1 , iz+3 , it , 4 , 1 ) = .true.
                  mask( ix+2 , iy+2 , iz+0 , it , 4 , 1 ) = .true.
                  mask( ix+2 , iy+3 , iz+1 , it , 4 , 1 ) = .true.
!     plane 4
                  mask( ix+3 , iy+0 , iz+3 , it , 4 , 1 ) = .true.
                  mask( ix+3 , iy+1 , iz+0 , it , 4 , 1 ) = .true.
                  mask( ix+3 , iy+2 , iz+1 , it , 4 , 1 ) = .true.
                  mask( ix+3 , iy+3 , iz+2 , it , 4 , 1 ) = .true.

               end do
            end do
         end do
      end do

!
!  Now shift this mask 1 step in t direction and 2 in x, y, and z directions
!    and OR it to the first
!
      mask(:,:,:,:,4,1) = mask(:,:,:,:,4,1) .or.
     &         cshift(
     &           cshift(
     &             cshift(
     &               cshift(mask(:,:,:,:,4,1), dim=4, shift=1),
     &             dim=1, shift=2),
     &           dim=2, shift=2),
     &         dim=3, shift=2)


!     mask 2, 3, 4

      do imask = 2, ncube
        mask(:,:,:,:,4,imask) = cshift( mask(:,:,:,:,4,imask-1), dim=1, shift=1 )
      end do

      return

      end subroutine MaskUltraFast

      END MODULE GS_MASKULTRAFAST
