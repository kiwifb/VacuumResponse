c     subroutine to calculate values for L shaped wilson loops WL, 
c     with "that" being the time direction
c     Adapted from wilsonloop.f and the original shapeTloop.f by F. Bissey Dec. 2003
c
c     It is also supplemented with the facilities to change the plane of computation:
c     xdir and ydir define the plane. 
c
c     This particular version takes links in the desired plane as input and delivers
c     4 L shaped Wilson Loops.
c     Last edited FB050324 
c
      MODULE L_LLOOPS

c     indexes corresponding to the steps in the xdir and ydir directions.
      integer                                                 :: ixq,iyq
      
      Contains

      subroutine l_loops(xdir,upxr,upxi,ydir,upyr,upyi,ur,ui,WL)

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

      end subroutine l_loops

c========================================================================================================

      END MODULE L_LLOOPS
