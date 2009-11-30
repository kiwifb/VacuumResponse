      MODULE GS_LATTICESIZE

      IMPLICIT NONE

!      integer, parameter :: nx=8,  ny=8,  nz=8,  nt=16
!      integer, parameter :: nx=6,  ny=6,  nz=6,  nt=12
      integer, parameter :: nx=12, ny=12, nz=12, nt=24
!      integer, parameter :: nx=16, ny=16, nz=16, nt=32
!      integer, parameter :: nx=24, ny=24, nz=24, nt=36
!      integer, parameter :: nx=20, ny=20, nz=20, nt=40
      integer,parameter  :: mMask=4
      integer,parameter  :: nIndex = nx*ny*nz*nt/mMask

      END MODULE GS_LATTICESIZE
