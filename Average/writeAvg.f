c  
c     ===============================================================
c     subroutine to write correlated action and topological charge 
c     from averaging program in unformated and formated file at a 
c     reasonnable speed on more than one node.
c
c     F. Bissey 18. Oct. 2004
c
      MODULE L_WRITEAVG

      CONTAINS

      subroutine unf_writeavg(handle,TopAct,Wavg,avgTopAct)

      USE L_baryonParam

      integer                                                 :: handle

      double precision,dimension(0:ny-1,0:nz-1,0:nt-1)        :: TopAct
!HPF$ PROCESSORS linProcs(NUMBER_OF_PROCESSORS())
!HPF$ TEMPLATE form(nt*NUMBER_OF_PROCESSORS())
!HPF$ ALIGN TopAct(*,*,:) WITH form(1:nt)
!HPF$ DISTRIBUTE form(BLOCK) onto linProcs
      double precision                                        :: Wavg
      double precision                                        :: avgTopAct
c
      write(handle) TopAct(:,:,:)
      write(handle) Wavg, avgTopAct

      return

      end subroutine unf_writeavg

      subroutine fmt_writeavg(handle,TopAct)

      USE L_baryonParam

      integer                                                 :: handle

      double precision,dimension(0:ny-1,0:nz-1,0:nt-1)        :: TopAct
!HPF$ PROCESSORS linProcs(NUMBER_OF_PROCESSORS())
!HPF$ TEMPLATE form(nt*NUMBER_OF_PROCESSORS())
!HPF$ ALIGN TopAct(*,*,:) WITH form(1:nt)
!HPF$ DISTRIBUTE form(BLOCK) onto linProcs

      write(handle,'(f15.8)') TopAct(:,:,:)

      return

      end subroutine fmt_writeavg

      END MODULE L_WRITEAVG
