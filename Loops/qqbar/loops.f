c
c     subroutine to calculate valules for Wilson loops, W, for 'Y' shape
c     v.3.11 AK040129
c
      MODULE L_LOOPS

      CONTAINS

c========================================================================================================

      subroutine qq_loops(xdir,ur,ui,W)
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
      integer                                                           :: xdir
      double precision,dimension(nx,ny,nz,nt,mu,nc,nc)                  :: ur,ui
!HPF$ DISTRIBUTE ur(*,*,BLOCK,BLOCK,*,*,*)
!HPF$ DISTRIBUTE ui(*,*,BLOCK,BLOCK,*,*,*)
c
      double precision,dimension(nx,ny,nz,nt,0:nloop,nL(that))          :: W
!HPF$ DISTRIBUTE W(*,*,BLOCK,BLOCK,*,*)
c
c
c     local variables
c
      integer                                                           :: it,iy
      integer                                                           :: ic,jc,kc
c
c     bottom and time links (built up iteratively)
c
      double precision,dimension(nx,ny,nz,nt,nc,nc)                     :: bxr,bxi
!HPF$ DISTRIBUTE bxr(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE bxi(*,*,BLOCK,BLOCK,*,*)
      double precision,dimension(nx,ny,nz,nt,nc,nc)                     :: tlr,tli
!HPF$ DISTRIBUTE tlr(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE tli(*,*,BLOCK,BLOCK,*,*)
c
c     shifted time links
c
      double precision,dimension(nx,ny,nz,nt,nc,nc)                     :: utr,uti
!HPF$ DISTRIBUTE utr(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE uti(*,*,BLOCK,BLOCK,*,*)
c
c      top links (bottom shifted)
c
      double precision,dimension(nx,ny,nz,nt,nc,nc)                     :: uxr,uxi
!HPF$ DISTRIBUTE uxr(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE uxi(*,*,BLOCK,BLOCK,*,*)
c
c     temporary arrays
c
      double precision,dimension(nx,ny,nz,nt,nc,nc)                     :: tmpr,tmpi,stmpr,stmpi
!HPF$ DISTRIBUTE tmpr(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE tmpi(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE stmpr(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE stmpi(*,*,BLOCK,BLOCK,*,*)
      double precision,dimension(nx,ny,nz,nt)                           :: Res
!HPF$ DISTRIBUTE Res(*,*,BLOCK,BLOCK)
c
c     ----------------
c     execution begins
c     ----------------
c
c     ---------------------------------
c     start of Wilson loop calculations
c     ---------------------------------
c
c     it do loop
c
      do it = 1, nL(that)

         call product(ur,ui,that,tlr,tli,it)

         do iy = 1, nloop

            call product(ur,ui,xdir,bxr,bxi,iy)

            utr = cshift(tlr, dim = xdir , shift = iy)
            uti = cshift(tli, dim = xdir , shift = iy)

            uxr = cshift(bxr, dim = that , shift = it)
            uxi = cshift(bxi, dim = that , shift = it)

            tmpr = 0.d0
            tmpi = 0.d0
            stmpr = 0.d0
            stmpi = 0.d0

            do ic = 1, nc
               do jc = 1, nc
                  do kc = 1, nc

                     tmpr(:,:,:,:,ic,jc) = tmpr(:,:,:,:,ic,jc) +
     &                  bxr(:,:,:,:,ic,kc)*utr(:,:,:,:,kc,jc) - bxi(:,:,:,:,ic,kc)*uti(:,:,:,:,kc,jc)

                     tmpi(:,:,:,:,ic,jc) = tmpi(:,:,:,:,ic,jc) +
     &                  bxr(:,:,:,:,ic,kc)*uti(:,:,:,:,kc,jc) + bxi(:,:,:,:,ic,kc)*utr(:,:,:,:,kc,jc)

                     stmpr(:,:,:,:,ic,jc) = stmpr(:,:,:,:,ic,jc) +
     &                  tlr(:,:,:,:,ic,kc)*uxr(:,:,:,:,kc,jc) - tli(:,:,:,:,ic,kc)*uxi(:,:,:,:,kc,jc)

                     stmpi(:,:,:,:,ic,jc) = stmpi(:,:,:,:,ic,jc) +
     &                  tlr(:,:,:,:,ic,kc)*uxi(:,:,:,:,kc,jc) + tli(:,:,:,:,ic,kc)*uxr(:,:,:,:,kc,jc)

                  end do
               end do
            end do

            res = 0.d0

            do ic = 1, nc
               do jc = 1, nc

                  res(:,:,:,:) = res(:,:,:,:) +
     &               tmpr(:,:,:,:,ic,jc)*stmpr(:,:,:,:,ic,jc) + tmpi(:,:,:,:,ic,jc)*stmpi(:,:,:,:,ic,jc)

               end do
            end do

            W(:,:,:,:,iy,it) = res(:,:,:,:)

         end do

      end do

      return

      end subroutine qq_loops

      END MODULE L_LOOPS
