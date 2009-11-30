!
!     ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!     subroutine to write to disk the link variables.
!
!     Author: Derek B. Leinweber
!
      MODULE GS_WRITELINKS

      CONTAINS

      subroutine WriteLinks(filename,ur,ui,nfig,beta,nxf,nyf,nzf,ntf,lastPlaq,plaqbarAvg,uzero)

      USE GS_LATTICESIZE
      USE GS_GETACTION
      USE GS_FIXSU3

      implicit none

!     global variables

      integer,parameter                                       :: nc=3
      integer,parameter                                       :: mu=4
      double precision,dimension(nx,ny,nz,nt,mu,nc,nc)        :: ur,ui
!HPF$ DISTRIBUTE ur(*,*,BLOCK,BLOCK,*,*,*)
!HPF$ DISTRIBUTE ui(*,*,BLOCK,BLOCK,*,*,*)
      character(len=*)                                        :: filename
      integer                                                 :: nfig,nxf,nyf,nzf,ntf
      double precision                                        :: beta,lastPlaq,plaqbarAvg,uzero

!     Local variables
!
!  The following will ensure the links ur0 and ui0 sit on a single processor
!
      double precision,dimension(nx,ny,nz,nt,mu,nc,nc)        :: ur0,ui0

!HPF$ PROCESSORS linProcs(NUMBER_OF_PROCESSORS())
!HPF$ TEMPLATE form(nt*NUMBER_OF_PROCESSORS())
!HPF$ ALIGN ur0(*,*,*,:,*,*,*) WITH form(1:nt)
!HPF$ ALIGN ui0(*,*,*,:,*,*,*) WITH form(1:nt)
!HPF$ DISTRIBUTE form(BLOCK) onto linProcs

      double precision,dimension(nx,ny,nz,nt)                 :: action
!HPF$ DISTRIBUTE action(*,*,BLOCK,BLOCK)

      integer                                                 :: ic
!
!  Timer Support
!
      INTEGER start_count, end_count, count_rate
      REAL    elapsed_time
!
!  Useage
!
!  CALL SYSTEM_CLOCK(start_count, count_rate)
!  CALL SYSTEM_CLOCK(end_count)
!  elapsed_time = REAL(end_count - start_count) / REAL(count_rate)
!

!     Execution begins

      CALL SYSTEM_CLOCK(start_count, count_rate)

      nfig = 0
      call fixsu3(ur,ui)

      ur0 = ur
      ui0 = ui

!     Update lastPlaq here for error checking: Set itype parameter to zero to get the plaquette

      call GetAction(ur,ui,action,lastPlaq,0,uzero)

      open(101,file=filename,form='unformatted',status='new',action='write')

      write(101) nfig,beta,nxf,nyf,nzf,ntf

      do ic = 1, nc-1
        write(101) ur0(:,:,:,:,:,ic,:)
        write(101) ui0(:,:,:,:,:,ic,:)
      end do

      write(101) lastPlaq, plaqbarAvg, uzero
      close(101)

      CALL SYSTEM_CLOCK(end_count)
      elapsed_time = REAL(end_count - start_count) / REAL(count_rate)
!
      if (elapsed_time .ne. 0.0) then
        write(*,'(/,a,f7.3,a,/)') 'The elapsed time is',elapsed_time,
     &    ' seconds for writing the links.'
      end if

      return

      end subroutine WriteLinks

      END MODULE GS_WRITELINKS
