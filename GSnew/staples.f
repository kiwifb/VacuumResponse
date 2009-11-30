!
!     cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!     the subroutine staples combines the wilson staples(squares plaquettes) and
!     the improved staples(the rectangles plaquettes) for one passed in xhat
!     direction. It returns the variables stapler and staplei.
!
!     Authors: Frederi! D.R. Bonnet & Derek B. Leinweber
!     Date: September 1998.
!
      MODULE GS_STAPLES

      CONTAINS

      subroutine staples(ur,ui,stapler,staplei,xhat,local,itype,uzero)

      USE GS_LATTICESIZE
      USE GS_SQUARES
      USE GS_RECTANGLES

      implicit none

!     global variables

      integer,parameter                                       :: nc=3
      integer,parameter                                       :: mu=4

      integer                                                 :: xhat
      integer                                                 :: itype
      double precision                                        :: uzero

      logical                                                 :: local

      double precision,dimension(nx,ny,nz,nt,nc,nc)           :: stapler,staplei
!HPF$ DISTRIBUTE stapler(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE staplei(*,*,BLOCK,BLOCK,*,*)
      double precision,dimension(nx,ny,nz,nt,mu,nc,nc)        :: ur,ui
!HPF$ DISTRIBUTE ur(*,*,BLOCK,BLOCK,*,*,*)
!HPF$ DISTRIBUTE ui(*,*,BLOCK,BLOCK,*,*,*)

!     local variables

      double precision,dimension(nx,ny,nz,nt,nc,nc)           :: rectr,recti        !rectangular staples
!HPF$ DISTRIBUTE rectr(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE recti(*,*,BLOCK,BLOCK,*,*)

!     start of the executions commands

      call squares(ur,ui,stapler,staplei,xhat,local)

      if(itype==1) then
         call rectangles(ur,ui,rectr,recti,xhat,local)
         stapler = ( 5.0d0 / 3.0d0 ) * stapler - ( 1.0d0 / ( 12.0d0 * uzero**2 ) ) * rectr
         staplei = ( 5.0d0 / 3.0d0 ) * staplei - ( 1.0d0 / ( 12.0d0 * uzero**2 ) ) * recti
      end if

      return

      end subroutine staples

      END MODULE GS_STAPLES
