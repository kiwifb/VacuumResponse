c
c  A program to merge old Y loops going from 0-7 with data going from 8-9
c  Adapted from AverageYall
c
      program merge89
c
      USE L_baryonParam
c
      IMPLICIT NONE
c
      integer,parameter                                       :: offmax=4
      integer,parameter                                       :: nYloop=9
      integer                                                 :: start, finish, icon, Correlate
      character(len=132)                                      :: configName, suffixName
      character(len=132)                                      :: thisConfig
      character(len=12)                                       :: fstr1,fstr2,wstr1,wstr2
      character(len=30)                                       :: dialog
c
      double precision,dimension(nL(that),offmax)             :: avgAction, avgTopChg
!HPF$ DISTRIBUTE avgAction(*,*)
!HPF$ DISTRIBUTE avgTopChg(*,*)
c
      double precision,dimension(0:ny-1,0:nz-1,0:nt-1,0:nYloop,nL(that),offmax)   :: actionC
!HPF$ DISTRIBUTE actionC(*,BLOCK,BLOCK,*,*,*)
      double precision,dimension(0:ny-1,0:nz-1,0:nt-1,0:nYloop,nL(that),offmax)   :: topChgC
!HPF$ DISTRIBUTE topChgC(*,BLOCK,BLOCK,*,*,*)
      double precision,dimension(0:nYloop,nL(that))           :: Wavg
!HPF$ DISTRIBUTE Wavg(*,*)
c
c  Correlations
c
      write(*,*)
      write(*,'(a)')'What would you like to correlate?'
      write(*,'(a)')'      1: Action and topological charge.'
      write(*,'(a)')'      2: Electric and Magnetic fields.'
      read (*,*) Correlate
c
      write(*,'(a)') 'What is the base file name?'
      read (*,'(a80)') configName
c
      write(*,'(a)') 'Please provide a range of configuration numbers. e.g. 0 100'
      read (*,*) start, finish
c
c  The suffix is of the form ".s 4" for example.
c
      write(*,'(a)') 'What is the file name suffix?'
      read (*,'(a80)') suffixName
c'
c Initialiaze file strings and dialogue strings
c
      if ( Correlate == 1) then
c
         fstr1='.action'
         fstr2='.topChg'
         wstr1='Action'
         wstr2='TopChg'
         dialog='action and topological '
c
      else
c
         fstr1='.ele'
         fstr2='.mag'
         wstr1='Electric'
         wstr2='Magnetic'
         dialog='electric and magnetic field '
c
      endif
c
      do icon = start, finish
c
c  Create file name
c
         write(thisConfig,fmt='(a,i3.3,a)') trim(configName),icon,trim(suffixName)
c         
         write(*,*) trim(thisConfig)
c
c  -------------------------------------
c  read in results for one configuration
c  -------------------------------------
c
         write(*,'(/,3a,/)') 'Reading ',dialog,'correlations.'
c
c  .Yp-xy
c
         open (11,file=thisConfig(1:len_trim(thisConfig))//'.Yp-xy'//fstr1(1:len_trim(fstr1))//'.correl.unf',
     &        status='old',form='unformatted')
         open (12,file=thisConfig(1:len_trim(thisConfig))//'.Yp-xy'//fstr2(1:len_trim(fstr2))//'.correl.unf',
     &        status='old',form='unformatted')
         open (13,file=thisConfig(1:len_trim(thisConfig))//'.Yp-xy89'//fstr1(1:len_trim(fstr1))//'.correl.unf',
     &        status='old',form='unformatted')
         open (14,file=thisConfig(1:len_trim(thisConfig))//'.Yp-xy89'//fstr2(1:len_trim(fstr2))//'.correl.unf',
     &        status='old',form='unformatted')
c
c  ----------
c  reading 0-7
c  ----------
c
         read(11) actionC(:,:,:,0:7,:,:)
         read(11) Wavg(0:7,:), avgAction(:,:)
         read(12) topChgC(:,:,:,0:7,:,:)
         read(12) Wavg(0:7,:), avgTopChg(:,:)
c
c  ----------
c  reading 8-9
c  ----------
c
         read(13) actionC(:,:,:,8:9,:,:)
         read(13) Wavg(8:9,:), avgAction(:,:)
         read(14) topChgC(:,:,:,8:9,:,:)
         read(14) Wavg(8:9,:), avgTopChg(:,:)
c
         close(11)
         close(12)
         close(13)
         close(14)
c
c  ---------------
c  writing 0-9
c  ---------------
c
         open (11,file=thisConfig(1:len_trim(thisConfig))//'.Yp-xy'//fstr1(1:len_trim(fstr1))//'.correl.unf',
     &        status='old',form='unformatted')
         open (12,file=thisConfig(1:len_trim(thisConfig))//'.Yp-xy'//fstr2(1:len_trim(fstr2))//'.correl.unf',
     &        status='old',form='unformatted')

         write(11) actionC(:,:,:,:,:,:)
         write(11) Wavg(:,:), avgAction(:,:)
         write(12) topChgC(:,:,:,:,:,:)
         write(12) Wavg(:,:), avgTopChg(:,:)
c
         close(11)
         close(12)
c
c  .Ym-xy
c
         open (11,file=thisConfig(1:len_trim(thisConfig))//'.Ym-xy'//fstr1(1:len_trim(fstr1))//'.correl.unf',
     &        status='old',form='unformatted')
         open (12,file=thisConfig(1:len_trim(thisConfig))//'.Ym-xy'//fstr2(1:len_trim(fstr2))//'.correl.unf',
     &        status='old',form='unformatted')
         open (13,file=thisConfig(1:len_trim(thisConfig))//'.Ym-xy89'//fstr1(1:len_trim(fstr1))//'.correl.unf',
     &        status='old',form='unformatted')
         open (14,file=thisConfig(1:len_trim(thisConfig))//'.Ym-xy89'//fstr2(1:len_trim(fstr2))//'.correl.unf',
     &        status='old',form='unformatted')
c
c  ----------
c  reading 0-7
c  ----------
c
         read(11) actionC(:,:,:,0:7,:,:)
         read(11) Wavg(0:7,:), avgAction(:,:)
         read(12) topChgC(:,:,:,0:7,:,:)
         read(12) Wavg(0:7,:), avgTopChg(:,:)
c
c  ----------
c  reading 8-9
c  ----------
c
         read(13) actionC(:,:,:,8:9,:,:)
         read(13) Wavg(8:9,:), avgAction(:,:)
         read(14) topChgC(:,:,:,8:9,:,:)
         read(14) Wavg(8:9,:), avgTopChg(:,:)
c
         close(11)
         close(12)
         close(13)
         close(14)
c
c  ---------------
c  writing 0-9
c  ---------------
c
         open (11,file=thisConfig(1:len_trim(thisConfig))//'.Ym-xy'//fstr1(1:len_trim(fstr1))//'.correl.unf',
     &        status='old',form='unformatted')
         open (12,file=thisConfig(1:len_trim(thisConfig))//'.Ym-xy'//fstr2(1:len_trim(fstr2))//'.correl.unf',
     &        status='old',form='unformatted')

         write(11) actionC(:,:,:,:,:,:)
         write(11) Wavg(:,:), avgAction(:,:)
         write(12) topChgC(:,:,:,:,:,:)
         write(12) Wavg(:,:), avgTopChg(:,:)
c
         close(11)
         close(12)
c
c  .Yp-xz 
c
         open (11,file=thisConfig(1:len_trim(thisConfig))//'.Yp-xz'//fstr1(1:len_trim(fstr1))//'.correl.unf',
     &        status='old',form='unformatted')
         open (12,file=thisConfig(1:len_trim(thisConfig))//'.Yp-xz'//fstr2(1:len_trim(fstr2))//'.correl.unf',
     &        status='old',form='unformatted')
         open (13,file=thisConfig(1:len_trim(thisConfig))//'.Yp-xz89'//fstr1(1:len_trim(fstr1))//'.correl.unf',
     &        status='old',form='unformatted')
         open (14,file=thisConfig(1:len_trim(thisConfig))//'.Yp-xz89'//fstr2(1:len_trim(fstr2))//'.correl.unf',
     &        status='old',form='unformatted')
c
c  ----------
c  reading 0-7
c  ----------
c
         read(11) actionC(:,:,:,0:7,:,:)
         read(11) Wavg(0:7,:), avgAction(:,:)
         read(12) topChgC(:,:,:,0:7,:,:)
         read(12) Wavg(0:7,:), avgTopChg(:,:)
c
c  ----------
c  reading 8-9
c  ----------
c
         read(13) actionC(:,:,:,8:9,:,:)
         read(13) Wavg(8:9,:), avgAction(:,:)
         read(14) topChgC(:,:,:,8:9,:,:)
         read(14) Wavg(8:9,:), avgTopChg(:,:)
c
         close(11)
         close(12)
         close(13)
         close(14)
c
c  ---------------
c  writing 0-9
c  ---------------
c
         open (11,file=thisConfig(1:len_trim(thisConfig))//'.Yp-xz'//fstr1(1:len_trim(fstr1))//'.correl.unf',
     &        status='old',form='unformatted')
         open (12,file=thisConfig(1:len_trim(thisConfig))//'.Yp-xz'//fstr2(1:len_trim(fstr2))//'.correl.unf',
     &        status='old',form='unformatted')

         write(11) actionC(:,:,:,:,:,:)
         write(11) Wavg(:,:), avgAction(:,:)
         write(12) topChgC(:,:,:,:,:,:)
         write(12) Wavg(:,:), avgTopChg(:,:)
c
         close(11)
         close(12)
c
c  .Ym-xz
c
         open (11,file=thisConfig(1:len_trim(thisConfig))//'.Ym-xz'//fstr1(1:len_trim(fstr1))//'.correl.unf',
     &        status='old',form='unformatted')
         open (12,file=thisConfig(1:len_trim(thisConfig))//'.Ym-xz'//fstr2(1:len_trim(fstr2))//'.correl.unf',
     &        status='old',form='unformatted')
         open (13,file=thisConfig(1:len_trim(thisConfig))//'.Ym-xz89'//fstr1(1:len_trim(fstr1))//'.correl.unf',
     &        status='old',form='unformatted')
         open (14,file=thisConfig(1:len_trim(thisConfig))//'.Ym-xz89'//fstr2(1:len_trim(fstr2))//'.correl.unf',
     &        status='old',form='unformatted')
c
c  ----------
c  reading 0-7
c  ----------
c
         read(11) actionC(:,:,:,0:7,:,:)
         read(11) Wavg(0:7,:), avgAction(:,:)
         read(12) topChgC(:,:,:,0:7,:,:)
         read(12) Wavg(0:7,:), avgTopChg(:,:)
c
c  ----------
c  reading 8-9
c  ----------
c
         read(13) actionC(:,:,:,8:9,:,:)
         read(13) Wavg(8:9,:), avgAction(:,:)
         read(14) topChgC(:,:,:,8:9,:,:)
         read(14) Wavg(8:9,:), avgTopChg(:,:)
c
         close(11)
         close(12)
         close(13)
         close(14)
c
c  ---------------
c  writing 0-9
c  ---------------
c
         open (11,file=thisConfig(1:len_trim(thisConfig))//'.Ym-xz'//fstr1(1:len_trim(fstr1))//'.correl.unf',
     &        status='old',form='unformatted')
         open (12,file=thisConfig(1:len_trim(thisConfig))//'.Ym-xz'//fstr2(1:len_trim(fstr2))//'.correl.unf',
     &        status='old',form='unformatted')

         write(11) actionC(:,:,:,:,:,:)
         write(11) Wavg(:,:), avgAction(:,:)
         write(12) topChgC(:,:,:,:,:,:)
         write(12) Wavg(:,:), avgTopChg(:,:)
c
         close(11)
         close(12)

      end do                    ! end icon loop

      end program merge89
