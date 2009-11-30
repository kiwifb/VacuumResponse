!
!     cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!     this subroutine reads in the link variables.
!
!     Author: Derek B. Leinweber
!
      MODULE GS_READLINKSRISKY

      CONTAINS

      subroutine ReadLinks(filename,ur,ui,nfig,beta,nxf,nyf,nzf,ntf,lastPlaq,plaqbarAvg,uzero)

      USE GS_LATTICESIZE
      USE GS_FIXSU3
      USE GS_GETACTION

      implicit none

!     global variables

      integer,parameter                                       :: nc=3
      integer,parameter                                       :: mu=4

      character(len=*)                                        :: filename
      integer                                                 :: nfig,nxf,nyf,nzf,ntf
      double precision                                        :: beta,lastPlaq,plaqbarAvg,uzero

      double precision,dimension(nx,ny,nz,nt,mu,nc,nc)        :: ur,ui
!HPF$ DISTRIBUTE ur(*,*,BLOCK,BLOCK,*,*,*)
!HPF$ DISTRIBUTE ui(*,*,BLOCK,BLOCK,*,*,*)

!     Variables to support GetAction

      double precision,dimension(nx,ny,nz,nt)                 :: action
!HPF$ DISTRIBUTE action(*,*,BLOCK,BLOCK)
      double precision                                        :: plaqbar


!     Local variables
!
!  The following will ensure the links ur0 and ui0 sit on two processors
!
      double precision,dimension(nx,ny,nz,nt,mu,nc,nc)        :: ur0,ui0

!HPF$ PROCESSORS linProcs(NUMBER_OF_PROCESSORS())
!HPF$ TEMPLATE form(nt*NUMBER_OF_PROCESSORS()/2)
!HPF$ ALIGN ur0(*,*,*,:,*,*,*) WITH form(1:nt)
!HPF$ ALIGN ui0(*,*,*,:,*,*,*) WITH form(1:nt)
!HPF$ DISTRIBUTE form(BLOCK) onto linProcs

      integer                                                 ::  ic

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

      write(*,'(/,a,a)') 'Opening file ',filename(1:len_trim(filename))
      open(101,file=filename,form='unformatted',status='old',action='read')

      write(*,'(a)') 'Reading preliminary data.'
      read (101) nfig,beta,nxf,nyf,nzf,ntf
      write(*,'(a,i4,a,f5.3)') 'Configuration number ',nfig,' at beta = ',beta
      write(*,'(a,4(i2,a))')   'Lattice Size ',nxf,'x',nyf,'x',nzf,'x',ntf

      write(*,'(a)') 'Reading the links.'
      do ic = 1, nc-1
        read (101) ur0(:,:,:,:,:,ic,:)
        read (101) ui0(:,:,:,:,:,ic,:)
      end do

      ur = ur0
      ui = ui0

      write(*,'(a)') 'Reading the final data.'
      read (101) lastPlaq, plaqbarAvg, uzero
      close(101)

      CALL SYSTEM_CLOCK(end_count)
      elapsed_time = REAL(end_count - start_count) / REAL(count_rate)
!
      if (elapsed_time .ne. 0.0) then
        write(*,'(/,a,f7.3,a,/)') 'The elapsed time is',elapsed_time,
     &    ' seconds for reading the links.'
      end if

      write(*,'(a,f9.6,a,f9.6)') 'The last plaquette and closing uzero is  ',
     &                           lastPlaq,' and ',uzero
      write(*,'(a,f9.6)')        'The average of the last 20 plaquettes is ',
     &                           plaqbarAvg

      write(*,'(a)') 'Call fix SU(3).'
      call fixsu3(ur,ui)
!
!     Check average plaquette
!
      write(*,'(a)') 'Calculating the action.'
      call GetAction(ur,ui,action,plaqbar,0,uzero)
      write(*,'(a,f19.16)') 'Plaquette action is ', plaqbar
      write(*,'(a,f19.16)') 'Last Plaquette was  ', lastPlaq
!
      write(*,'(/,a)') 'Error checking suppressed in readLinksRisky'

      return

      end subroutine ReadLinks

      END MODULE GS_READLINKSRISKY
