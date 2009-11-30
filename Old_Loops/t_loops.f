c     subroutine to calculate values for T shaped wilson loops, WT*, 
c     with "that" being the time direction
c     Adapted from wilsonloop.f and the original shapeTloop.f by F. Bissey Dec. 2003
c     This particular version initially limited itself to 2 kinds of T shapes
c     2 kinds of L shapes. For statistics it has been decided to go back to 4 L shapes
c     as they are all fully exploitables.
c
c     It is also supplemented with the facilities to change the plane of computation:
c     xdir and ydir define the plane. 
c
      MODULE L_TLOOPS

c     indexes corresponding to the steps in the xdir and ydir directions.
      integer                                                 :: ixq,iyq

      Contains

      subroutine t_loops(xdir,upxr,upxi,ydir,upyr,upyi,ur,ui,WTp,WTm)

        USE L_baryonParam
        USE L_epsilonIndex
        USE L_product

        implicit none

        integer                                               :: xdir,ydir
        double precision,dimension(nx,ny,nz,nt,mu,nc,nc)      :: ur,ui        !links
!HPF$ DISTRIBUTE ur(*,*,BLOCK,BLOCK,*,*,*)
!HPF$ DISTRIBUTE ui(*,*,BLOCK,BLOCK,*,*,*)
c     links products in xdir and ydir directions
      double precision,dimension(nx,ny,nz,nt,nc,nc)           :: upxr,upxi,upyr,upyi
!HPF$ DISTRIBUTE upxr(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE upxi(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE upyr(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE upyi(*,*,BLOCK,BLOCK,*,*)
c     T-shapes are labelled p(plus) and m(minus):
c     p: ( ixq; iyq;-iyq)
c     m: (-ixq; iyq;-iyq)
      double precision,dimension(nx,ny,nz,nt,nL(that))        :: WTp,WTm ! T-shape wilson loops
!HPF$ DISTRIBUTE WTp(*,*,BLOCK,BLOCK,*)
!HPF$ DISTRIBUTE WTm(*,*,BLOCK,BLOCK,*)
c
      include 'LTfiles/Common_declarations.f'
c
c  Execution Begins
c
        do it = 1,nL(that)

           include 'LTfiles/Common_lt.f'
c
c==========================================================================
c==========================================================================
c==========================================================================
c
c  Construct the baryon for quarks at +iyq, -iyq and +ixq : WTp
c
           Res = 0.0d0

           do ic = 1,36

              Res(:,:,:,:) = Res(:,:,:,:)  + fper(ic) * (
     &               sypr(:,:,:,:,la(ic),lap(ic)) * synr(:,:,:,:,lb(ic),lbp(ic)) * sxpr(:,:,:,:,lc(ic),lcp(ic))
     &             - sypr(:,:,:,:,la(ic),lap(ic)) * syni(:,:,:,:,lb(ic),lbp(ic)) * sxpi(:,:,:,:,lc(ic),lcp(ic))
     &             - sypi(:,:,:,:,la(ic),lap(ic)) * synr(:,:,:,:,lb(ic),lbp(ic)) * sxpi(:,:,:,:,lc(ic),lcp(ic))
     &             - sypi(:,:,:,:,la(ic),lap(ic)) * syni(:,:,:,:,lb(ic),lbp(ic)) * sxpr(:,:,:,:,lc(ic),lcp(ic))
     &             )

           end do
c
           WTp(:,:,:,:,it) = Res(:,:,:,:) / 6.0d0
c
c  Construct the baryon for quarks at +iyq, -iyq and -ixq : WTm
c
           Res = 0.0d0

           do ic = 1,36

              Res(:,:,:,:) = Res(:,:,:,:)  + fper(ic) * (
     &               sypr(:,:,:,:,la(ic),lap(ic)) * synr(:,:,:,:,lb(ic),lbp(ic)) * sxnr(:,:,:,:,lc(ic),lcp(ic))
     &             - sypr(:,:,:,:,la(ic),lap(ic)) * syni(:,:,:,:,lb(ic),lbp(ic)) * sxni(:,:,:,:,lc(ic),lcp(ic))
     &             - sypi(:,:,:,:,la(ic),lap(ic)) * synr(:,:,:,:,lb(ic),lbp(ic)) * sxni(:,:,:,:,lc(ic),lcp(ic))
     &             - sypi(:,:,:,:,la(ic),lap(ic)) * syni(:,:,:,:,lb(ic),lbp(ic)) * sxnr(:,:,:,:,lc(ic),lcp(ic))
     &             )

           end do
c
           WTm(:,:,:,:,it) = Res(:,:,:,:)  / 6.0d0

        end do                  !end L and T loop

        return

      end subroutine t_loops

c========================================================================================================

      END MODULE L_TLOOPS
