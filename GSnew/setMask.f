!
!
!------------------------------------------------------------------------------------------
!
!  Subroutine SetMask sets the mask for a variety of actions and lattice sizes.
!
!  Author: Derek B. Leinweber  dleinweb@physics.adelaide.edu.au
!          14 Sept. 2000
!
!  For the standard Wilson action, call the routine with itype = 0.
!
!  For plaquette plus rectangle action, call the routine with itype = 1 or 2.
!    The routine will select the appropriate mask based on the lattice dimensions.
!
!  For 3, 4, or 5-loop improved actions, call the routine with itype = 3 or greater.
!
!  Upon return the mask will be set, and umask will be set to the number of masks
!    per link direction required to cover all links.
!
      MODULE GS_SETMASK

      CONTAINS

      subroutine SetMask(mask, umask, itype)

      USE GS_LATTICESIZE
      USE GS_MASKWILSON
      USE GS_MASKRECT
      USE GS_MASKULTRAFAST
      USE GS_MASKULTRA
!
      implicit none

      integer,parameter                       :: mu=4
      integer,parameter                       :: nmask=16
      integer                                 :: umask, itype
      logical,dimension(nx,ny,nz,nt,mu,nmask) :: mask
!HPF$ DISTRIBUTE mask(*,*,BLOCK,BLOCK,*,*)
!
      if (itype == 0) then
        call MaskWilson(mask)
        umask = 2
        write(*,'(/,a)')   'Selecting the 2x2 Standard Wilson Mask.'
        write(*,'(i2,a,/)') umask, ' masks per link direction are required.'
      else if (itype == 1 .or. itype == 2) then
!
!  Try for UltraFastMask first
!
        if      ( (nx/4)*4 == nx .and.
     &            (ny/4)*4 == ny .and.
     &            (nz/4)*4 == nz .and.
     &            (nt/4)*4 == nt       ) then
          call MaskUltraFast(mask)
          umask = 4
          write(*,'(/,a)')   'Selecting the 2x4 Fast Mask for Plaquette plus Rectangle.'
          write(*,'(i2,a,/)') umask, ' masks per link direction are required.'
        else if ( (nx/6)*6 == nx .and.
     &            (ny/6)*6 == ny .and.
     &            (nz/6)*6 == nz .and.
     &            (nt/6)*6 == nt       ) then
          call MaskRect(mask)
          umask = 6
          write(*,'(/,a)')   'Selecting the 2x3 Mask for Plaquette plus Rectangle.'
          write(*,'(i2,a,/)') umask, ' masks per link direction are required.'
        else
          write(*,'(/,a,/)') 'Unable to set the mask for these lattice dimensions.'
          mask = .false.
          umask = 0
        end if
!
      else if (itype >= 3) then
!
!  The lattice must be a multiple of 4
!
        if      ( (nx/4)*4 == nx .and.
     &            (ny/4)*4 == ny .and.
     &            (nz/4)*4 == nz .and.
     &            (nt/4)*4 == nt       ) then
!
!  If it is also a multiple of six we can do it with 12 masks
!
          if ( (nx/6)*6 == nx .and.
     &         (ny/6)*6 == ny .and.
     &         (nz/6)*6 == nz .and.
     &         (nt/6)*6 == nt       ) then
            call MaskUltra(mask,3)
            umask = 12
            write(*,'(/,a)')   'Selecting the 3x4 Mask for the O(a^4) Improved Action.'
            write(*,'(i2,a,/)') umask, ' masks per link direction are required.'
          else
            call MaskUltra(mask,4)
            umask = 16
            write(*,'(/,a)')   'Selecting the 4x4 Mask for the O(a^4) Improved Action.'
            write(*,'(i2,a,/)') umask, ' masks per link direction are required.'
          end if
        else
          write(*,'(/,a,/)') 'Unable to set the mask for these lattice dimensions.'
          mask = .false.
          umask = 0
        end if
      else
        write(*,'(/,a,/)') 'Invalid action type.'
        mask = .false.
        umask = 0
      end if
!
      return

      end subroutine SetMask

      END MODULE GS_SETMASK
