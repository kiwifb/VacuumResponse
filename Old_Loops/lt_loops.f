c     subroutine to calculate values for T and L shaped wilson loops, WT* and WL*, 
c     with "that" being the time direction
c     Adapted from wilsonloop.f and the original shapeTloop.f by F. Bissey Dec. 2003
c     This particular version initially limited itself to 2 kinds of T shapes and
c     2 kinds of L shapes. For statistics it has been decided to go back to 4 L shapes
c     as they are all fully exploitables.
c
c     It is also supplemented with the facilities to change the plane of computation:
c     xdir and ydir define the plane. 
c
c     20/07/2004 Merged lt-loops and and lt-loops-redux adding plane configurations
c     capabilities to both. lt_loops compute 2 T and 4 L shapes, lt_lopps_all 
c     computes 4 T and 4 L shapes. Also their common parts have been separated and 
c     put in LTfiles in a similar way to y-loops. F. Bissey.
c
      MODULE L_LTLOOPS

c     indexes corresponding to the steps in the xdir and ydir directions.
      integer                                                 :: ixq,iyq

      Contains

      subroutine lt_loops(xdir,upxr,upxi,ydir,upyr,upyi,ur,ui,WTp,WTm,WL)

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
c     L-shapes are numbered 1 to 4 :
c     1: ( ixq, iyq)
c     2: ( ixq,-iyq)
c     3: (-ixq,-iyq)
c     4: (-ixq, iyq)
      double precision,dimension(nx,ny,nz,nt,nL(that),4)      :: WL ! L-shape wilson loops
!HPF$ DISTRIBUTE WL(*,*,BLOCK,BLOCK,*,*)
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
c
c  Construct the baryon for quarks at +iyq, 0 and +ixq : WL(1)
c
           Res = 0.0d0

           do ic = 1,36

              Res(:,:,:,:) = Res(:,:,:,:)  + fper(ic) * (
     &               sypr(:,:,:,:,la(ic),lap(ic)) * uptr(:,:,:,:,lb(ic),lbp(ic)) * sxpr(:,:,:,:,lc(ic),lcp(ic))
     &             - sypr(:,:,:,:,la(ic),lap(ic)) * upti(:,:,:,:,lb(ic),lbp(ic)) * sxpi(:,:,:,:,lc(ic),lcp(ic))
     &             - sypi(:,:,:,:,la(ic),lap(ic)) * uptr(:,:,:,:,lb(ic),lbp(ic)) * sxpi(:,:,:,:,lc(ic),lcp(ic))
     &             - sypi(:,:,:,:,la(ic),lap(ic)) * upti(:,:,:,:,lb(ic),lbp(ic)) * sxpr(:,:,:,:,lc(ic),lcp(ic))
     &             )

           end do
c
           WL(:,:,:,:,it,1) = Res(:,:,:,:) / 6.0d0
c
c  Construct the baryon for quarks at -iyq, 0 and ixq : WL(2)
c
           Res = 0.0d0

           do ic = 1,36

              Res(:,:,:,:) = Res(:,:,:,:)  + fper(ic) * (
     &               synr(:,:,:,:,la(ic),lap(ic)) * uptr(:,:,:,:,lb(ic),lbp(ic)) * sxpr(:,:,:,:,lc(ic),lcp(ic))
     &             - synr(:,:,:,:,la(ic),lap(ic)) * upti(:,:,:,:,lb(ic),lbp(ic)) * sxpi(:,:,:,:,lc(ic),lcp(ic))
     &             - syni(:,:,:,:,la(ic),lap(ic)) * uptr(:,:,:,:,lb(ic),lbp(ic)) * sxpi(:,:,:,:,lc(ic),lcp(ic))
     &             - syni(:,:,:,:,la(ic),lap(ic)) * upti(:,:,:,:,lb(ic),lbp(ic)) * sxpr(:,:,:,:,lc(ic),lcp(ic))
     &             )

           end do
c
           WL(:,:,:,:,it,2) = Res(:,:,:,:)  / 6.0d0
c
c  Construct the baryon for quarks at -iyq, 0 and -ixq : WL(3)
c
           Res = 0.0d0

           do ic = 1,36

              Res(:,:,:,:) = Res(:,:,:,:)  + fper(ic) * (
     &               synr(:,:,:,:,la(ic),lap(ic)) * uptr(:,:,:,:,lb(ic),lbp(ic)) * sxnr(:,:,:,:,lc(ic),lcp(ic))
     &             - synr(:,:,:,:,la(ic),lap(ic)) * upti(:,:,:,:,lb(ic),lbp(ic)) * sxni(:,:,:,:,lc(ic),lcp(ic))
     &             - syni(:,:,:,:,la(ic),lap(ic)) * uptr(:,:,:,:,lb(ic),lbp(ic)) * sxni(:,:,:,:,lc(ic),lcp(ic))
     &             - syni(:,:,:,:,la(ic),lap(ic)) * upti(:,:,:,:,lb(ic),lbp(ic)) * sxnr(:,:,:,:,lc(ic),lcp(ic))
     &             )

           end do
c
           WL(:,:,:,:,it,3) = Res(:,:,:,:)  / 6.0d0
c
c  Construct the baryon for quarks at +iyq, 0 and -ixq : WL(4)
c
           Res = 0.0d0

           do ic = 1,36

              Res(:,:,:,:) = Res(:,:,:,:)  + fper(ic) * (
     &               sypr(:,:,:,:,la(ic),lap(ic)) * uptr(:,:,:,:,lb(ic),lbp(ic)) * sxnr(:,:,:,:,lc(ic),lcp(ic))
     &             - sypr(:,:,:,:,la(ic),lap(ic)) * upti(:,:,:,:,lb(ic),lbp(ic)) * sxni(:,:,:,:,lc(ic),lcp(ic))
     &             - sypi(:,:,:,:,la(ic),lap(ic)) * uptr(:,:,:,:,lb(ic),lbp(ic)) * sxni(:,:,:,:,lc(ic),lcp(ic))
     &             - sypi(:,:,:,:,la(ic),lap(ic)) * upti(:,:,:,:,lb(ic),lbp(ic)) * sxnr(:,:,:,:,lc(ic),lcp(ic))
     &             )

           end do
c
           WL(:,:,:,:,it,4) = Res(:,:,:,:)  / 6.0d0

        end do                  !end L and T loop

        return

      end subroutine lt_loops

c========================================================================================================

      subroutine lt_loops_all(xdir,upxr,upxi,ydir,upyr,upyi,ur,ui,WLpp,WLpm,WLmp,WLmm,WTpp,WTpm,WTmp,WTmm)

        USE L_baryonParam
        USE L_epsilonIndex
        USE L_product

        implicit none

        integer                                                 :: xdir,ydir
        double precision,dimension(nx,ny,nz,nt,mu,nc,nc)        :: ur,ui        !links
!HPF$ DISTRIBUTE ur(*,*,BLOCK,BLOCK,*,*,*)
!HPF$ DISTRIBUTE ui(*,*,BLOCK,BLOCK,*,*,*)
        double precision,dimension(nx,ny,nz,nt,nc,nc)           :: upxr,upxi,upyr,upyi
!HPF$ DISTRIBUTE upxr(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE upxi(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE upyr(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE upyi(*,*,BLOCK,BLOCK,*,*)
        double precision,dimension(nx,ny,nz,nt,nL(that))        :: WLpp,WLpm,WLmp,WLmm ! L-shape wilson loops
!HPF$ DISTRIBUTE WLpp(*,*,BLOCK,BLOCK,*)
!HPF$ DISTRIBUTE WLpm(*,*,BLOCK,BLOCK,*)
!HPF$ DISTRIBUTE WLmp(*,*,BLOCK,BLOCK,*)
!HPF$ DISTRIBUTE WLmm(*,*,BLOCK,BLOCK,*)
        double precision,dimension(nx,ny,nz,nt,nL(that))        :: WTpp,WTpm,WTmp,WTmm ! T-shape wilson loops
!HPF$ DISTRIBUTE WTpp(*,*,BLOCK,BLOCK,*)
!HPF$ DISTRIBUTE WTpm(*,*,BLOCK,BLOCK,*)
!HPF$ DISTRIBUTE WTmp(*,*,BLOCK,BLOCK,*)
!HPF$ DISTRIBUTE WTmm(*,*,BLOCK,BLOCK,*)
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
c  Construct the baryon for quarks at +iyq, -iyq and +ixq : WTpp
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
           WTpp(:,:,:,:,it) = Res(:,:,:,:) / 6.0d0
c
c  Construct the baryon for quarks at +iyq, -iyq and -ixq : WTmm
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
           WTmm(:,:,:,:,it) = Res(:,:,:,:)  / 6.0d0
c
c  Construct the baryon for quarks at +ixq, -ixq and +iyq : WTmp
c
           Res = 0.0d0

           do ic = 1,36

              Res(:,:,:,:) = Res(:,:,:,:)  + fper(ic) * (
     &               sxpr(:,:,:,:,la(ic),lap(ic)) * sxnr(:,:,:,:,lb(ic),lbp(ic)) * sypr(:,:,:,:,lc(ic),lcp(ic))
     &             - sxpr(:,:,:,:,la(ic),lap(ic)) * sxni(:,:,:,:,lb(ic),lbp(ic)) * sypi(:,:,:,:,lc(ic),lcp(ic))
     &             - sxpi(:,:,:,:,la(ic),lap(ic)) * sxnr(:,:,:,:,lb(ic),lbp(ic)) * sypi(:,:,:,:,lc(ic),lcp(ic))
     &             - sxpi(:,:,:,:,la(ic),lap(ic)) * sxni(:,:,:,:,lb(ic),lbp(ic)) * sypr(:,:,:,:,lc(ic),lcp(ic))
     &             )

           end do
c
           WTmp(:,:,:,:,it) = Res(:,:,:,:) / 6.0d0
c
c  Construct the baryon for quarks at +ixq, -ixq and -iyq : WTpm
c
           Res = 0.0d0

           do ic = 1,36

              Res(:,:,:,:) = Res(:,:,:,:)  + fper(ic) * (
     &               sxpr(:,:,:,:,la(ic),lap(ic)) * sxnr(:,:,:,:,lb(ic),lbp(ic)) * synr(:,:,:,:,lc(ic),lcp(ic))
     &             - sxpr(:,:,:,:,la(ic),lap(ic)) * sxni(:,:,:,:,lb(ic),lbp(ic)) * syni(:,:,:,:,lc(ic),lcp(ic))
     &             - sxpi(:,:,:,:,la(ic),lap(ic)) * sxnr(:,:,:,:,lb(ic),lbp(ic)) * syni(:,:,:,:,lc(ic),lcp(ic))
     &             - sxpi(:,:,:,:,la(ic),lap(ic)) * sxni(:,:,:,:,lb(ic),lbp(ic)) * synr(:,:,:,:,lc(ic),lcp(ic))
     &             )

           end do
c
           WTpm(:,:,:,:,it) = Res(:,:,:,:)  / 6.0d0
c
c  Construct the baryon for quarks at +iyq, 0 and +ixq : WLpp
c
           Res = 0.0d0

           do ic = 1,36

              Res(:,:,:,:) = Res(:,:,:,:)  + fper(ic) * (
     &               sypr(:,:,:,:,la(ic),lap(ic)) * uptr(:,:,:,:,lb(ic),lbp(ic)) * sxpr(:,:,:,:,lc(ic),lcp(ic))
     &             - sypr(:,:,:,:,la(ic),lap(ic)) * upti(:,:,:,:,lb(ic),lbp(ic)) * sxpi(:,:,:,:,lc(ic),lcp(ic))
     &             - sypi(:,:,:,:,la(ic),lap(ic)) * uptr(:,:,:,:,lb(ic),lbp(ic)) * sxpi(:,:,:,:,lc(ic),lcp(ic))
     &             - sypi(:,:,:,:,la(ic),lap(ic)) * upti(:,:,:,:,lb(ic),lbp(ic)) * sxpr(:,:,:,:,lc(ic),lcp(ic))
     &             )

           end do
c
           WLpp(:,:,:,:,it) = Res(:,:,:,:) / 6.0d0
c
c  Construct the baryon for quarks at +iyq, 0 and -ixq : WLmp
c
           Res = 0.0d0

           do ic = 1,36

              Res(:,:,:,:) = Res(:,:,:,:)  + fper(ic) * (
     &               sypr(:,:,:,:,la(ic),lap(ic)) * uptr(:,:,:,:,lb(ic),lbp(ic)) * sxnr(:,:,:,:,lc(ic),lcp(ic))
     &             - sypr(:,:,:,:,la(ic),lap(ic)) * upti(:,:,:,:,lb(ic),lbp(ic)) * sxni(:,:,:,:,lc(ic),lcp(ic))
     &             - sypi(:,:,:,:,la(ic),lap(ic)) * uptr(:,:,:,:,lb(ic),lbp(ic)) * sxni(:,:,:,:,lc(ic),lcp(ic))
     &             - sypi(:,:,:,:,la(ic),lap(ic)) * upti(:,:,:,:,lb(ic),lbp(ic)) * sxnr(:,:,:,:,lc(ic),lcp(ic))
     &             )

           end do
c
           WLmp(:,:,:,:,it) = Res(:,:,:,:)  / 6.0d0
c
c  Construct the baryon for quarks at -iyq, 0 and ixq : WLpm
c
           Res = 0.0d0

           do ic = 1,36

              Res(:,:,:,:) = Res(:,:,:,:)  + fper(ic) * (
     &               synr(:,:,:,:,la(ic),lap(ic)) * uptr(:,:,:,:,lb(ic),lbp(ic)) * sxpr(:,:,:,:,lc(ic),lcp(ic))
     &             - synr(:,:,:,:,la(ic),lap(ic)) * upti(:,:,:,:,lb(ic),lbp(ic)) * sxpi(:,:,:,:,lc(ic),lcp(ic))
     &             - syni(:,:,:,:,la(ic),lap(ic)) * uptr(:,:,:,:,lb(ic),lbp(ic)) * sxpi(:,:,:,:,lc(ic),lcp(ic))
     &             - syni(:,:,:,:,la(ic),lap(ic)) * upti(:,:,:,:,lb(ic),lbp(ic)) * sxpr(:,:,:,:,lc(ic),lcp(ic))
     &             )

           end do
c
           WLpm(:,:,:,:,it) = Res(:,:,:,:)  / 6.0d0
c
c  Construct the baryon for quarks at -iyq, 0 and -ixq : WLmm
c
           Res = 0.0d0

           do ic = 1,36

              Res(:,:,:,:) = Res(:,:,:,:)  + fper(ic) * (
     &               synr(:,:,:,:,la(ic),lap(ic)) * uptr(:,:,:,:,lb(ic),lbp(ic)) * sxnr(:,:,:,:,lc(ic),lcp(ic))
     &             - synr(:,:,:,:,la(ic),lap(ic)) * upti(:,:,:,:,lb(ic),lbp(ic)) * sxni(:,:,:,:,lc(ic),lcp(ic))
     &             - syni(:,:,:,:,la(ic),lap(ic)) * uptr(:,:,:,:,lb(ic),lbp(ic)) * sxni(:,:,:,:,lc(ic),lcp(ic))
     &             - syni(:,:,:,:,la(ic),lap(ic)) * upti(:,:,:,:,lb(ic),lbp(ic)) * sxnr(:,:,:,:,lc(ic),lcp(ic))
     &             )

           end do
c
           WLmm(:,:,:,:,it) = Res(:,:,:,:)  / 6.0d0
c

        end do                  !end L and T loop

        return

      end subroutine lt_loops_all

      END MODULE L_LTLOOPS
