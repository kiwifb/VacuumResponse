c
c     subroutines to calculate sets of fixed Wilson loops, W, for 'L' shape
c
      MODULE L_LOOPS

      CONTAINS

c========================================================================================================

      subroutine ly_loops(xdir,ydir,ur,ui,uptr,upti,it,W)
c
c     modules
c
      USE L_baryonParam
      USE L_epsilonIndex
      USE L_product
c
c     no implicit typing
c
      implicit none
c
c     global variables
c
#include"loopsize.f"

      integer                                                           :: xdir,ydir
      double precision,dimension(nx,ny,nz,nt,mu,nc,nc)                  :: ur,ui
!HPF$ DISTRIBUTE ur(*,*,BLOCK,BLOCK,*,*,*)
!HPF$ DISTRIBUTE ui(*,*,BLOCK,BLOCK,*,*,*)
      double precision,dimension(nx,ny,nz,nt,nc,nc)                     :: uptr,upti
!HPF$ DISTRIBUTE uptr(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE upti(*,*,BLOCK,BLOCK,*,*)
c
c            4|1
c            -+-
c            3|2
c
      double precision,dimension(nx,ny,nz,nt,0:nloop,4)                 :: W
!HPF$ DISTRIBUTE W(*,*,BLOCK,BLOCK,*,*)
c
c     local variables
c
      integer                                                           :: it,illoop
      integer                                                           :: ic,jc,kc
c
c     bottom and top links (built up iteratively)
c
      double precision,dimension(nx,ny,nz,nt,nc,nc)                     :: bpxlr,bpxli
!HPF$ DISTRIBUTE bpxlr(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE bpxli(*,*,BLOCK,BLOCK,*,*)
      double precision,dimension(nx,ny,nz,nt,nc,nc)                     :: bpylr,bpyli
!HPF$ DISTRIBUTE bpylr(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE bpyli(*,*,BLOCK,BLOCK,*,*)
      double precision,dimension(nx,ny,nz,nt,nc,nc)                     :: bmxlr,bmxli
!HPF$ DISTRIBUTE bmxlr(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE bmxli(*,*,BLOCK,BLOCK,*,*)
      double precision,dimension(nx,ny,nz,nt,nc,nc)                     :: bmylr,bmyli
!HPF$ DISTRIBUTE bmylr(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE bmyli(*,*,BLOCK,BLOCK,*,*)
      double precision,dimension(nx,ny,nz,nt,nc,nc)                     :: tlr,tli
!HPF$ DISTRIBUTE tlr(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE tli(*,*,BLOCK,BLOCK,*,*)
c
c     shifted time links
c
      double precision,dimension(nx,ny,nz,nt,nc,nc)                     :: ustr,usti
!HPF$ DISTRIBUTE ustr(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE usti(*,*,BLOCK,BLOCK,*,*)
c
c     staples/temporary arrays
c
      double precision,dimension(nx,ny,nz,nt,nc,nc)                     :: spxlr,spxli
!HPF$ DISTRIBUTE spxlr(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE spxli(*,*,BLOCK,BLOCK,*,*)
      double precision,dimension(nx,ny,nz,nt,nc,nc)                     :: spylr,spyli
!HPF$ DISTRIBUTE spylr(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE spyli(*,*,BLOCK,BLOCK,*,*)
      double precision,dimension(nx,ny,nz,nt,nc,nc)                     :: smxlr,smxli
!HPF$ DISTRIBUTE smxlr(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE smxli(*,*,BLOCK,BLOCK,*,*)
      double precision,dimension(nx,ny,nz,nt,nc,nc)                     :: smylr,smyli
!HPF$ DISTRIBUTE smylr(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE smyli(*,*,BLOCK,BLOCK,*,*)
c
c     temporary arrays
c
      double precision,dimension(nx,ny,nz,nt,nc,nc)                     :: tmpr,tmpi
!HPF$ DISTRIBUTE tmpr(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE tmpi(*,*,BLOCK,BLOCK,*,*)
      double precision,dimension(nx,ny,nz,nt)                           :: Res
!HPF$ DISTRIBUTE Res(*,*,BLOCK,BLOCK)
c
c     ---------------------------------
c     start of Wilson loop calculations
c     ---------------------------------
      do illoop = 1, nloop

         call product(ur,ui,xdir,bpxlr,bpxli,illoop)
         call product(ur,ui,ydir,bpylr,bpyli,illoop)

         bmxlr = cshift(bpxlr,dim=xdir,shift=-illoop)
         bmxli = cshift(bpxli,dim=xdir,shift=-illoop)
         bmylr = cshift(bpylr,dim=ydir,shift=-illoop)
         bmyli = cshift(bpyli,dim=ydir,shift=-illoop)
c
c     Make +xdir staple
c
         ustr(:,:,:,:,:,:) = cshift(uptr(:,:,:,:,:,:),dim=xdir,shift=illoop)
         usti(:,:,:,:,:,:) = cshift(upti(:,:,:,:,:,:),dim=xdir,shift=illoop)

         tmpr = 0.0d0
         tmpi = 0.0d0

         do ic = 1,nc
            do jc = 1,nc
               do kc = 1,nc

                  tmpr(:,:,:,:,ic,jc) = tmpr(:,:,:,:,ic,jc)
     &                 + bpxlr(:,:,:,:,ic,kc) * ustr(:,:,:,:,kc,jc) - bpxli(:,:,:,:,ic,kc) * usti(:,:,:,:,kc,jc)

                  tmpi(:,:,:,:,ic,jc) = tmpi(:,:,:,:,ic,jc)
     &                 + bpxlr(:,:,:,:,ic,kc) * usti(:,:,:,:,kc,jc) + bpxli(:,:,:,:,ic,kc) * ustr(:,:,:,:,kc,jc)

               end do
            end do
         end do

         spxlr = 0.0d0
         spxli = 0.0d0
         tlr = cshift(bpxlr,dim=that,shift=it)
         tli = cshift(bpxli,dim=that,shift=it)

         do ic = 1,nc
            do jc = 1,nc
               do kc = 1,nc

                  spxlr(:,:,:,:,ic,jc) = spxlr(:,:,:,:,ic,jc)
     &                 + tmpr(:,:,:,:,ic,kc) * tlr(:,:,:,:,jc,kc) + tmpi(:,:,:,:,ic,kc) * tli(:,:,:,:,jc,kc)

                  spxli(:,:,:,:,ic,jc) = spxli(:,:,:,:,ic,jc)
     &                 - tmpr(:,:,:,:,ic,kc) * tli(:,:,:,:,jc,kc) + tmpi(:,:,:,:,ic,kc) * tlr(:,:,:,:,jc,kc)

               end do
            end do
         end do
c
c     Make +ydir staple
c
         ustr(:,:,:,:,:,:) = cshift(uptr(:,:,:,:,:,:),dim=ydir,shift=illoop)
         usti(:,:,:,:,:,:) = cshift(upti(:,:,:,:,:,:),dim=ydir,shift=illoop)

         tmpr = 0.0d0
         tmpi = 0.0d0

         do ic = 1,nc
            do jc = 1,nc
               do kc = 1,nc

                  tmpr(:,:,:,:,ic,jc) = tmpr(:,:,:,:,ic,jc)
     &                 + bpylr(:,:,:,:,ic,kc) * ustr(:,:,:,:,kc,jc) - bpyli(:,:,:,:,ic,kc) * usti(:,:,:,:,kc,jc)

                  tmpi(:,:,:,:,ic,jc) = tmpi(:,:,:,:,ic,jc)
     &                 + bpylr(:,:,:,:,ic,kc) * usti(:,:,:,:,kc,jc) + bpyli(:,:,:,:,ic,kc) * ustr(:,:,:,:,kc,jc)

               end do
            end do
         end do

         spylr = 0.0d0
         spyli = 0.0d0
         tlr = cshift(bpylr,dim=that,shift=it)
         tli = cshift(bpyli,dim=that,shift=it)

         do ic = 1,nc
            do jc = 1,nc
               do kc = 1,nc

                  spylr(:,:,:,:,ic,jc) = spylr(:,:,:,:,ic,jc)
     &                 + tmpr(:,:,:,:,ic,kc) * tlr(:,:,:,:,jc,kc) + tmpi(:,:,:,:,ic,kc) * tli(:,:,:,:,jc,kc)

                  spyli(:,:,:,:,ic,jc) = spyli(:,:,:,:,ic,jc)
     &                 - tmpr(:,:,:,:,ic,kc) * tli(:,:,:,:,jc,kc) + tmpi(:,:,:,:,ic,kc) * tlr(:,:,:,:,jc,kc)

               end do
            end do
         end do
c
c     Make -xdir staple
c
         ustr(:,:,:,:,:,:) = cshift(uptr(:,:,:,:,:,:),dim=xdir,shift=-illoop)
         usti(:,:,:,:,:,:) = cshift(upti(:,:,:,:,:,:),dim=xdir,shift=-illoop)

         tmpr = 0.0d0
         tmpi = 0.0d0

         do ic = 1,nc
            do jc = 1,nc
               do kc = 1,nc

                  tmpr(:,:,:,:,ic,jc) = tmpr(:,:,:,:,ic,jc)
     &                 + bmxlr(:,:,:,:,kc,ic) * ustr(:,:,:,:,kc,jc) + bmxli(:,:,:,:,kc,ic) * usti(:,:,:,:,kc,jc)

                  tmpi(:,:,:,:,ic,jc) = tmpi(:,:,:,:,ic,jc)
     &                 + bmxlr(:,:,:,:,kc,ic) * usti(:,:,:,:,kc,jc) - bmxli(:,:,:,:,kc,ic) * ustr(:,:,:,:,kc,jc)

               end do
            end do
         end do

         smxlr = 0.0d0
         smxli = 0.0d0
         tlr = cshift(bmxlr,dim=that,shift=it)
         tli = cshift(bmxli,dim=that,shift=it)

         do ic = 1,nc
            do jc = 1,nc
               do kc = 1,nc

                  smxlr(:,:,:,:,ic,jc) = smxlr(:,:,:,:,ic,jc)
     &                 + tmpr(:,:,:,:,ic,kc) * tlr(:,:,:,:,kc,jc) - tmpi(:,:,:,:,ic,kc) * tli(:,:,:,:,kc,jc)

                  smxli(:,:,:,:,ic,jc) = smxli(:,:,:,:,ic,jc)
     &                 + tmpr(:,:,:,:,ic,kc) * tli(:,:,:,:,kc,jc) + tmpi(:,:,:,:,ic,kc) * tlr(:,:,:,:,kc,jc)

               end do
            end do
         end do
c
c     Make -ydir staple
c
         ustr(:,:,:,:,:,:) = cshift(uptr(:,:,:,:,:,:),dim=ydir,shift=-illoop)
         usti(:,:,:,:,:,:) = cshift(upti(:,:,:,:,:,:),dim=ydir,shift=-illoop)

         tmpr = 0.0d0
         tmpi = 0.0d0

         do ic = 1,nc
            do jc = 1,nc
               do kc = 1,nc

                  tmpr(:,:,:,:,ic,jc) = tmpr(:,:,:,:,ic,jc)
     &                 + bmylr(:,:,:,:,kc,ic) * ustr(:,:,:,:,kc,jc) + bmyli(:,:,:,:,kc,ic) * usti(:,:,:,:,kc,jc)

                  tmpi(:,:,:,:,ic,jc) = tmpi(:,:,:,:,ic,jc)
     &                 + bmylr(:,:,:,:,kc,ic) * usti(:,:,:,:,kc,jc) - bmyli(:,:,:,:,kc,ic) * ustr(:,:,:,:,kc,jc)

               end do
            end do
         end do

         smylr = 0.0d0
         smyli = 0.0d0
         tlr = cshift(bmylr,dim=that,shift=it)
         tli = cshift(bmyli,dim=that,shift=it)

         do ic = 1,nc
            do jc = 1,nc
               do kc = 1,nc

                  smylr(:,:,:,:,ic,jc) = smylr(:,:,:,:,ic,jc)
     &                 + tmpr(:,:,:,:,ic,kc) * tlr(:,:,:,:,kc,jc) - tmpi(:,:,:,:,ic,kc) * tli(:,:,:,:,kc,jc)

                  smyli(:,:,:,:,ic,jc) = smyli(:,:,:,:,ic,jc)
     &                 + tmpr(:,:,:,:,ic,kc) * tli(:,:,:,:,kc,jc) + tmpi(:,:,:,:,ic,kc) * tlr(:,:,:,:,kc,jc)

               end do
            end do
         end do
c
c     Make Wilson loop 1 |_
c
         Res = 0.0d0
         do ic = 1,36

            Res(:,:,:,:) = Res(:,:,:,:)  + fper(ic) * (
     &           spxlr(:,:,:,:,la(ic),lap(ic)) * spylr(:,:,:,:,lb(ic),lbp(ic)) * uptr(:,:,:,:,lc(ic),lcp(ic))
     &         - spxlr(:,:,:,:,la(ic),lap(ic)) * spyli(:,:,:,:,lb(ic),lbp(ic)) * upti(:,:,:,:,lc(ic),lcp(ic))
     &         - spxli(:,:,:,:,la(ic),lap(ic)) * spylr(:,:,:,:,lb(ic),lbp(ic)) * upti(:,:,:,:,lc(ic),lcp(ic))
     &         - spxli(:,:,:,:,la(ic),lap(ic)) * spyli(:,:,:,:,lb(ic),lbp(ic)) * uptr(:,:,:,:,lc(ic),lcp(ic))
     &           )

         end do

         W(:,:,:,:,illoop,1) = Res(:,:,:,:)
c                         _
c     Make Wilson loop 2 |
c
         Res = 0.0d0
         do ic = 1,36

            Res(:,:,:,:) = Res(:,:,:,:)  + fper(ic) * (
     &           spxlr(:,:,:,:,la(ic),lap(ic)) * smylr(:,:,:,:,lb(ic),lbp(ic)) * uptr(:,:,:,:,lc(ic),lcp(ic))
     &         - spxlr(:,:,:,:,la(ic),lap(ic)) * smyli(:,:,:,:,lb(ic),lbp(ic)) * upti(:,:,:,:,lc(ic),lcp(ic))
     &         - spxli(:,:,:,:,la(ic),lap(ic)) * smylr(:,:,:,:,lb(ic),lbp(ic)) * upti(:,:,:,:,lc(ic),lcp(ic))
     &         - spxli(:,:,:,:,la(ic),lap(ic)) * smyli(:,:,:,:,lb(ic),lbp(ic)) * uptr(:,:,:,:,lc(ic),lcp(ic))
     &           )

         end do

         W(:,:,:,:,illoop,2) = Res(:,:,:,:)
c                         _
c     Make Wilson loop 3   |
c
         Res = 0.0d0
         do ic = 1,36

            Res(:,:,:,:) = Res(:,:,:,:)  + fper(ic) * (
     &           smxlr(:,:,:,:,la(ic),lap(ic)) * smylr(:,:,:,:,lb(ic),lbp(ic)) * uptr(:,:,:,:,lc(ic),lcp(ic))
     &         - smxlr(:,:,:,:,la(ic),lap(ic)) * smyli(:,:,:,:,lb(ic),lbp(ic)) * upti(:,:,:,:,lc(ic),lcp(ic))
     &         - smxli(:,:,:,:,la(ic),lap(ic)) * smylr(:,:,:,:,lb(ic),lbp(ic)) * upti(:,:,:,:,lc(ic),lcp(ic))
     &         - smxli(:,:,:,:,la(ic),lap(ic)) * smyli(:,:,:,:,lb(ic),lbp(ic)) * uptr(:,:,:,:,lc(ic),lcp(ic))
     &           )

         end do

         W(:,:,:,:,illoop,3) = Res(:,:,:,:)
c
c     Make Wilson loop 4  _|
c
         Res = 0.0d0
         do ic = 1,36

            Res(:,:,:,:) = Res(:,:,:,:)  + fper(ic) * (
     &           smxlr(:,:,:,:,la(ic),lap(ic)) * spylr(:,:,:,:,lb(ic),lbp(ic)) * uptr(:,:,:,:,lc(ic),lcp(ic))
     &         - smxlr(:,:,:,:,la(ic),lap(ic)) * spyli(:,:,:,:,lb(ic),lbp(ic)) * upti(:,:,:,:,lc(ic),lcp(ic))
     &         - smxli(:,:,:,:,la(ic),lap(ic)) * spylr(:,:,:,:,lb(ic),lbp(ic)) * upti(:,:,:,:,lc(ic),lcp(ic))
     &         - smxli(:,:,:,:,la(ic),lap(ic)) * spyli(:,:,:,:,lb(ic),lbp(ic)) * uptr(:,:,:,:,lc(ic),lcp(ic))
     &           )

         end do

         W(:,:,:,:,illoop,4) = Res(:,:,:,:)

      end do                    ! end illoop

      return

      end subroutine ly_loops

      END MODULE L_LOOPS
