!
!
!------------------------------------------------------------------------------------------
!
!
!     Subroutine MaskRect establishes masks for the plaquette plus rectangle action
!       on lattice dimensions divisible by 6 but not 4.
!
!     Authors: Frederi! Bonnet           fbonnet@physics.adelaide.edu.au
!              Sundance Bilson-Thompson  sbilson@physics.adelaide.edu.au
!              Derek B. Leinweber        dleinweb@physics.adelaide.edu.au
!              14 Sept. 2000
!
      MODULE GS_MASKRECT

      CONTAINS

      subroutine MaskRect(mask)

      USE GS_LATTICESIZE

      implicit none
!     global variables

      integer,parameter                       :: mu=4
      integer,parameter                       :: nmask=16
      logical,dimension(nx,ny,nz,nt,mu,nmask) :: mask
!HPF$ DISTRIBUTE mask(*,*,BLOCK,BLOCK,*,*)

!     local variables

      integer,parameter                       :: ncube=3
      integer,parameter                       :: nlink=2
      integer,parameter                       :: nspace=3
      integer                                 :: ix,iy,iz,it
      integer                                 :: imask
      integer                                 :: idir, ilink

!
!------------------------------------------------------------------------------------------
!
      mask = .false.
!
!-----------------------------------------------------------------------------------------
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
!     plane 2
                  mask( ix , iy+0 , iz+1 , it+1 , 1 , 1 ) = .true.
                  mask( ix , iy+1 , iz+2 , it+1 , 1 , 1 ) = .true.
                  mask( ix , iy+2 , iz+0 , it+1 , 1 , 1 ) = .true.
!     plane 3
                  mask( ix , iy+0 , iz+2 , it+2 , 1 , 1 ) = .true.
                  mask( ix , iy+1 , iz+0 , it+2 , 1 , 1 ) = .true.
                  mask( ix , iy+2 , iz+1 , it+2 , 1 , 1 ) = .true.

               end do
            end do
         end do
      end do

!     mask 1 for the y direction

      do ix=1,nx,ncube
         do iy=1,ny,nlink
            do iz=1,nz,ncube
               do it=1,nt,ncube
!     plane 1
                  mask( ix+0 , iy , iz+0 , it+0 , 2 , 1 ) = .true.
                  mask( ix+1 , iy , iz+1 , it+0 , 2 , 1 ) = .true.
                  mask( ix+2 , iy , iz+2 , it+0 , 2 , 1 ) = .true.
!     plane 2
                  mask( ix+0 , iy , iz+1 , it+1 , 2 , 1 ) = .true.
                  mask( ix+1 , iy , iz+2 , it+1 , 2 , 1 ) = .true.
                  mask( ix+2 , iy , iz+0 , it+1 , 2 , 1 ) = .true.
!     plane 3
                  mask( ix+0 , iy , iz+2 , it+2 , 2 , 1 ) = .true.
                  mask( ix+1 , iy , iz+0 , it+2 , 2 , 1 ) = .true.
                  mask( ix+2 , iy , iz+1 , it+2 , 2 , 1 ) = .true.

               end do
            end do
         end do
      end do

!     mask 1 for the z direction

      do ix=1,nx,ncube
         do iy=1,ny,ncube
            do iz=1,nz,nlink
               do it=1,nt,ncube
!     plane 1
                  mask( ix+0 , iy+0 , iz , it+0 , 3 , 1 ) = .true.
                  mask( ix+1 , iy+1 , iz , it+0 , 3 , 1 ) = .true.
                  mask( ix+2 , iy+2 , iz , it+0 , 3 , 1 ) = .true.
!     plane 2
                  mask( ix+0 , iy+1 , iz , it+1 , 3 , 1 ) = .true.
                  mask( ix+1 , iy+2 , iz , it+1 , 3 , 1 ) = .true.
                  mask( ix+2 , iy+0 , iz , it+1 , 3 , 1 ) = .true.
!     plane 3
                  mask( ix+0 , iy+2 , iz , it+2 , 3 , 1 ) = .true.
                  mask( ix+1 , iy+0 , iz , it+2 , 3 , 1 ) = .true.
                  mask( ix+2 , iy+1 , iz , it+2 , 3 , 1 ) = .true.

               end do
            end do
         end do
      end do

!     mask 2, 3

      do idir = 1, nspace
        do imask = 2, ncube
          mask(:,:,:,:,idir,imask) = cshift( mask(:,:,:,:,idir,imask-1), dim=4, shift=1 )
        end do

!     mask 4 is mask 1 shifted by 1 in the x direction and so on

        do ilink = 1, nlink-1
          do imask = 1, ncube
            mask(:,:,:,:,idir,imask+ilink*ncube) =
     &          cshift( mask(:,:,:,:,idir,imask+(ilink-1)*ncube), dim=idir, shift=1 )
          end do
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
!     plane 2
                  mask( ix+1 , iy+0 , iz+1 , it , 4 , 1 ) = .true.
                  mask( ix+1 , iy+1 , iz+2 , it , 4 , 1 ) = .true.
                  mask( ix+1 , iy+2 , iz+0 , it , 4 , 1 ) = .true.
!     plane 3
                  mask( ix+2 , iy+0 , iz+2 , it , 4 , 1 ) = .true.
                  mask( ix+2 , iy+1 , iz+0 , it , 4 , 1 ) = .true.
                  mask( ix+2 , iy+2 , iz+1 , it , 4 , 1 ) = .true.

               end do
            end do
         end do
      end do

!     mask 2, 3

      do imask = 2, ncube
        mask(:,:,:,:,4,imask) = cshift( mask(:,:,:,:,4,imask-1), dim=1, shift=1 )
      end do

!     mask 4 is mask 1 shifted by 1 in the t direction and so on

      do ilink = 1, nlink-1
        do imask = 1, ncube
          mask(:,:,:,:,4,imask+ilink*ncube) =
     &        cshift( mask(:,:,:,:,4,imask+(ilink-1)*ncube), dim=4, shift=1 )
        end do
      end do

      return

      end subroutine MaskRect

      END MODULE GS_MASKRECT
