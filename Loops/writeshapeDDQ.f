c
c     ===============================================================
c     subroutine to write correlated action and topological charge
c     form Yshape in unformated file at a reasonnable speed on more
c     than one node.
c
c     F. Bissey 18 Aug. 2004
c
      MODULE L_WRITESHAPE

      CONTAINS

      subroutine writeshape(filename,TopAct,Wavg,avgTopAct,offmax)

      USE L_baryonParam

#include"loopsize.f"
      character(len=*)                                        :: filename
      integer                                                 :: offmax

      double precision,dimension(0:ny-1,0:nz-1,0:nt-1,1:ysep,1:dqsep,nL(that),offmax) :: TopAct
!HPF$ PROCESSORS linProcs(NUMBER_OF_PROCESSORS())
!HPF$ TEMPLATE form(nL(that)*NUMBER_OF_PROCESSORS())
!HPF$ ALIGN TopAct(*,*,*,*,:,:,*) WITH form(1:nL(that))
!HPF$ DISTRIBUTE form(BLOCK) onto linProcs
      double precision,dimension(1:ysep,1:dqsep,nL(that))      :: Wavg
!HPF$ DISTRIBUTE Wavg(*,*)
      double precision,dimension(nL(that),offmax)              :: avgTopAct
!HPF$ DISTRIBUTE avgTopAct(*,*)


      open(101,file=filename,form='unformatted',status='unknown',action='write')
c
      write(101) TopAct(:,:,:,:,:,:,:)
      write(101) Wavg(:,:,:), avgTopAct(:,:)

      close(101)

      return

      end subroutine writeshape

      END MODULE L_WRITESHAPE
