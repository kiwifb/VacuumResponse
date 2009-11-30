c  
c     ===============================================================
c     subroutine to write correlated action and topological charge 
c     form LTshape in unformated file at a reasonnable speed on more 
c     than one node.
c
c     F. Bissey 25 Aug. 2004
c
      MODULE L_WRITELTSHAPE

      CONTAINS

      subroutine writeLTshape(handle,TopAct,Wavg,avgTopAct,offmax)

      USE L_baryonParam

      integer                                                 :: handle
      integer                                                 :: offmax

      double precision,dimension(0:ny-1,0:nz-1,0:nt-1,nL(that),offmax) :: TopAct
!HPF$ PROCESSORS linProcs(NUMBER_OF_PROCESSORS())
!HPF$ TEMPLATE form(nL(that)*NUMBER_OF_PROCESSORS())
!HPF$ ALIGN TopAct(*,*,*,:,*) WITH form(1:nL(that))
!HPF$ DISTRIBUTE form(BLOCK) onto linProcs
      double precision,dimension(nL(that))                     :: Wavg
!HPF$ DISTRIBUTE Wavg(*)
      double precision,dimension(nL(that),offmax)              :: avgTopAct
!HPF$ DISTRIBUTE avgTopAct(*,*)


      write(handle) TopAct(:,:,:,:,:)
      write(handle) Wavg(:), avgTopAct(:,:)

      return

      end subroutine writeLTshape

      END MODULE L_WRITELTSHAPE 
