c
c  A program to test gauge invariance results.
c  Adapted from UniversalAverage
c
c  v.0.9  FB071121 Intial import
c
      program GaugeTest
c
      IMPLICIT NONE
c
      integer,parameter                                       :: nx=12, ny=12, nz=12, nt=24
      integer,parameter                                       :: nc=3         !sigma,colour
      integer,parameter                                       :: mu=4         !directions
      integer,parameter                                       :: xhat=4       !spatial direction for wilson loop
      integer,parameter                                       :: yhat=3       !spatial direction for wilson loop
      integer,parameter                                       :: zhat=2       !spatial direction for wilson loop
      integer,parameter                                       :: that=1       !time direction for wilson loop
      integer,parameter                                       :: Tmax=8
      integer,parameter                                       :: maxsize=12
      integer                                                 :: nloop
      integer,parameter                                       :: offmax=4
      integer                                                 :: Correlate
      integer                                                 :: offset,Tee
      character(len=80)                                       :: suffixName
      character(len=132)                                      :: thisConfig
      character(len=12)                                       :: fstr
      character(len=30)                                       :: dialog
c
      double precision                                        :: tmp
      double precision,dimension(Tmax,offmax)                 :: avgAction
!HPF$ DISTRIBUTE avgAction(BLOCK,BLOCK)
c
      double precision,dimension(0:ny-1,0:nz-1,0:nt-1,0:maxsize,Tmax,offmax)   :: actionC
!HPF$ DISTRIBUTE actionC(*,BLOCK,BLOCK,*,*,*)
      double precision,dimension(0:maxsize,Tmax)              :: Wavg
!HPF$ DISTRIBUTE Wavg(BLOCK,BLOCK)
      double precision,dimension(0:ny-1,0:nz-1,0:nt-1,0:maxsize,Tmax,offmax)   :: CAction_c
      double precision,dimension(0:ny-1,0:nz-1,0:nt-1,0:maxsize,Tmax,offmax)   :: CAction_rtgc
      integer                                                 :: ix,iy,iz,iiy
c
c  Begin execution
c
c  Basic input
c
      write(*,*)
      write(*,'(a)')'What would you like to test?'
      write(*,'(a)')'      1: Action'
      write(*,'(a)')'      2: topological charge.'
      write(*,'(a)')'      3: Electric field'
      write(*,'(a)')'      4: Magnetic field.'
      read (*,*) Correlate
c
c  The suffix is of the form ".s 4.TYavg4" for example.
c
      write(*,'(a)') 'What is the file name _long_ suffix?'
      read (*,'(a80)') suffixName
c
      write(*,'(a)') 'What is the number of loops?'
      read (*,*) nloop
c
c Initialiaze file strings and dialogue strings
c
      select case (Correlate)
      case(1)
         fstr='.action'
         dialog='action'
c
      case(2)
         fstr='.topChg'
         dialog='topological '
c
      case(3)
         fstr='.ele'
         dialog='electric field'
c
      case(4)
         fstr='.mag'
         dialog='magnetic field '
c
      end select
c
c  Create file name
c
         write(thisConfig,fmt='(a,a)') 'su3b460s12t24IMPc001',trim(suffixName)
c
         write(*,*) trim(thisConfig)
c
c  -------------------------------------
c  read in results for the c configuration
c  -------------------------------------
c
         write(*,'(3a)') 'Reading ',dialog,'correlations.'
c
         open (11,file=trim(thisConfig)//trim(fstr)//'.correl.unf',status='old',form='unformatted')
c
c
         read(11) actionC(:,:,:,0:nloop,:,:)
         read(11) Wavg(0:nloop,:), avgAction(:,:)
c
         close(11)
c
c  --------------
c  Computing the correlated action
c  --------------
c
         do Tee = 2,Tmax
            do offset = 1,offmax
c
               if ( Tee - 2*offset .lt. 0 ) then
                  cycle
               endif
c
               do iiy = 1, nloop
c
                  CAction_c(:,:,:,iiy,Tee,offset) =
     &                 actionC(:,:,:,iiy,Tee,offset) / (Wavg(iiy,Tee) * avgAction(Tee,offset))
c
               end do
            end do
         end do
c
c  Create file name
c
         write(thisConfig,fmt='(a,a)') 'su3b460s12t24IMPrgtc001',trim(suffixName)
c
         write(*,*) trim(thisConfig)
c
c  -------------------------------------
c  read in results for rtgc configuration
c  -------------------------------------
c
         write(*,'(3a)') 'Reading ',dialog,'correlations.'
c
         open (11,file=trim(thisConfig)//trim(fstr)//'.correl.unf',status='old',form='unformatted')
c
c
         read(11) actionC(:,:,:,0:nloop,:,:)
         read(11) Wavg(0:nloop,:), avgAction(:,:)
c
         close(11)
c
c  --------------
c  Computing the correlated action
c  --------------
c
         do Tee = 2,Tmax
            do offset = 1,offmax
c
               if ( Tee - 2*offset .lt. 0 ) then
                  cycle
               endif
c
               do iiy = 1, nloop
c
                  CAction_rtgc(:,:,:,iiy,Tee,offset) =
     &                 actionC(:,:,:,iiy,Tee,offset) / (Wavg(iiy,Tee) * avgAction(Tee,offset))
c
               end do
            end do
         end do
c
c
c  comparison results
c
         open (11,file='comp-'//trim(suffixName)//trim(fstr)//'.correl',status='unknown')
c
         do Tee = 2,Tmax
            do offset = 1,offmax
c
               if ( Tee - 2*offset .lt. 0 ) then
                  cycle
               endif
c
               do iiy = 1, nloop
c
                  write(11,'(a)')    '*********************'
                  write(11,'(a,3i3)') 'Quark separation scheme ', iiy, Tee,offset
                  write(11,'(a)')    '*********************'

                  do ix = 0, nt-1
                     do iy = 0, nz-1
                        do iz = 0, ny-1

                           tmp = abs(CAction_c(iz,iy,ix,iiy,Tee,offset)-CAction_rtgc(iz,iy,ix,iiy,Tee,offset))
                           if(tmp.ge.1e-7) then
                              write(11,'(f15.8,f15.8,f10.7)') CAction_c(iz,iy,ix,iiy,Tee,offset),
     &                           CAction_rtgc(iz,iy,ix,iiy,Tee,offset),tmp
                           end if

                        end do
                     end do
                  end do
c
               enddo               ! end iiy loop
c
            enddo                  ! end offset loop
         enddo                     ! end Tee loop
c
         close (11)
c
      end program GaugeTest
