c
c     subroutine to calculate valules for Wilson loops, W, for 'Quad quarks' shape
c     v.1 FB091228
c
      MODULE L_LOOPS

      CONTAINS

c========================================================================================================

      subroutine link_conjug(ur,ui)
c
c     modules
c
      USE GS_LATTICESIZE
      USE L_baryonParam
c
c     no implicit typing
c
      implicit none
c
c     global variables
c
      double precision,dimension(nx,ny,nz,nt,nc,nc)                  :: ur,ui
!HPF$ DISTRIBUTE ur(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE ui(*,*,BLOCK,BLOCK,*,*)
c
c     local variables
c
      double precision,dimension(nx,ny,nz,nt,nc,nc)                     :: tmpr,tmpi
!HPF$ DISTRIBUTE tmpr(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE tmpi(*,*,BLOCK,BLOCK,*,*)
      integer                                                           :: ic,jc

      do ic = 1, nc
         do jc = 1, nc

            tmpr(:,:,:,:,ic,jc) = ur(:,:,:,:,jc,ic)
            tmpi(:,:,:,:,ic,jc) =-ui(:,:,:,:,jc,ic)

         end do
      end do

      ur = tmpr
      ui = tmpi

      return

      end subroutine link_conjug

c========================================================================================================

      subroutine quad_loops(xdir,ydir,ur,ui,W)
c
c     modules
c
      USE GS_LATTICESIZE
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
c
      double precision,dimension(nx,ny,nz,nt,0:nloop,nL(that))          :: W
!HPF$ DISTRIBUTE W(*,*,BLOCK,BLOCK,*,*)
c
c
c     local variables
c
      integer                                                           :: it,ix,iy,iz
      integer                                                           :: lt,ie
      integer                                                           :: ic,jc,kc,mc
c
c     labelling of the staples and time links
c         a
c     1 _____ 2
c   d    | |     b
c     4 _|_|_ 3
c         c
c
c     time links (built up iteratively)
c     shifted time links
c
      double precision,dimension(nx,ny,nz,nt,nc,nc)                     :: t1r,t1i,t2r,t2i,t3r,t3i,t4r,t4i
!HPF$ DISTRIBUTE t1r(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE t1i(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE t2r(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE t2i(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE t3r(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE t3i(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE t4r(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE t4i(*,*,BLOCK,BLOCK,*,*)
c
c     bottom edge of the staples
c
      double precision,dimension(nx,ny,nz,nt,nc,nc)                     :: bx1r,bx1i,bx2r,bx2i,bx3r,bx3i,bx4r,bx4i
!HPF$ DISTRIBUTE bx1r(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE bx1i(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE bx2r(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE bx2i(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE bx3r(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE bx3i(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE bx4r(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE bx4i(*,*,BLOCK,BLOCK,*,*)
c
c     top edge of the staples
c
      double precision,dimension(nx,ny,nz,nt,nc,nc)                     :: ux1r,ux1i,ux2r,ux2i,ux3r,ux3i,ux4r,ux4i
!HPF$ DISTRIBUTE ux1r(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE ux1i(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE ux2r(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE ux2i(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE ux3r(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE ux3i(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE ux4r(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE ux4i(*,*,BLOCK,BLOCK,*,*)
c
c     staples
c
      double precision,dimension(nx,ny,nz,nt,nc,nc)                     :: s1r,s1i,s2r,s2i,s3r,s3i,s4r,s4i
!HPF$ DISTRIBUTE s1r(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE s1i(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE s2r(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE s2i(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE s3r(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE s3i(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE s4r(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE s4i(*,*,BLOCK,BLOCK,*,*)
c
c     Edges
c
      double precision,dimension(nx,ny,nz,nt,nc,nc)                     :: bar,bai,bbr,bbi,bcr,bci,bdr,bdi
!HPF$ DISTRIBUTE bar(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE bai(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE bbr(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE bbi(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE bcr(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE bci(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE bdr(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE bdi(*,*,BLOCK,BLOCK,*,*)
c
c     time shifted Edges
c
      double precision,dimension(nx,ny,nz,nt,nc,nc)                     :: uar,uai,ubr,ubi,ucr,uci,udr,udi
!HPF$ DISTRIBUTE uar(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE uai(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE ubr(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE ubi(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE ucr(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE uci(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE udr(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE udi(*,*,BLOCK,BLOCK,*,*)
c
c     "corners"
c
      double precision,dimension(nc,nc,nc,nc)              :: c1r,c1i,c2r,c2i
!HPF$ DISTRIBUTE c1r(BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE c1i(BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE c2r(BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE c2i(BLOCK,BLOCK,*,*)
c
c     temporary arrays
c
      double precision,dimension(nx,ny,nz,nt,nc,nc)                     :: tmp1r,tmp1i,tmp2r,tmp2i,tmp3r,tmp3i,tmp4r,tmp4i
!HPF$ DISTRIBUTE tmp1r(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE tmp1i(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE tmp2r(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE tmp2i(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE tmp3r(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE tmp3i(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE tmp4r(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE tmp4i(*,*,BLOCK,BLOCK,*,*)
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
c     staples edges
c
      call product(ur,ui,xdir,bx4r,bx4i,1)

      bx4r = cshift(bx4r, dim = xdir , shift =-1)
      bx4i = cshift(bx4i, dim = xdir , shift =-1)
      bx1r = cshift(bx4r, dim = ydir , shift = 2)
      bx1i = cshift(bx4i, dim = ydir , shift = 2)
      bx2r = cshift(bx1r, dim = xdir , shift = 2)
      bx2i = cshift(bx1i, dim = xdir , shift = 2)
      bx3r = cshift(bx4r, dim = xdir , shift = 2)
      bx3i = cshift(bx4i, dim = xdir , shift = 2)

      call product(ur,ui,xdir,bcr,bci,1)
      call product(ur,ui,xdir,bcr,bci,2)

      bar = cshift(bcr, dim = ydir , shift = 2)
      bai = cshift(bci, dim = ydir , shift = 2)

      call product(ur,ui,ydir,bdr,bdi,1)
      call product(ur,ui,ydir,bdr,bdi,2)

      bbr = cshift(bdr, dim = xdir , shift = 2)
      bbi = cshift(bdi, dim = xdir , shift = 2)
c
c     it do loop
c
      do lt = 1, nL(that)

         call product(ur,ui,that,t4r,t4i,lt)

         t1r = cshift(t4r, dim = ydir , shift = 2)
         t1i = cshift(t4i, dim = ydir , shift = 2)
         t2r = cshift(t1r, dim = xdir , shift = 2)
         t2i = cshift(t1i, dim = xdir , shift = 2)
         t3r = cshift(t4r, dim = xdir , shift = 2)
         t3i = cshift(t4i, dim = xdir , shift = 2)
c
c       Making the staples
c
         ux4r = cshift(bx4r, dim = that , shift =lt)
         ux4i = cshift(bx4i, dim = that , shift =lt)
         ux1r = cshift(bx1r, dim = that , shift =lt)
         ux1i = cshift(bx1i, dim = that , shift =lt)
         ux2r = cshift(bx2r, dim = that , shift =lt)
         ux2i = cshift(bx2i, dim = that , shift =lt)
         ux3r = cshift(bx3r, dim = that , shift =lt)
         ux3i = cshift(bx3i, dim = that , shift =lt)

         tmp1r = 0.d0
         tmp1i = 0.d0
         tmp2r = 0.d0
         tmp2i = 0.d0
         tmp3r = 0.d0
         tmp3i = 0.d0
         tmp4r = 0.d0
         tmp4i = 0.d0

         do ic = 1, nc
            do jc = 1, nc
               do kc = 1, nc

                  tmp4r(:,:,:,:,ic,jc) = tmp4r(:,:,:,:,ic,jc) +
     &               bx4r(:,:,:,:,kc,ic)*t4r(:,:,:,:,kc,jc) + bx4i(:,:,:,:,kc,ic)*t4i(:,:,:,:,kc,jc)

                  tmp4i(:,:,:,:,ic,jc) = tmp4i(:,:,:,:,ic,jc) +
     &               bx4r(:,:,:,:,kc,ic)*t4i(:,:,:,:,kc,jc) - bx4i(:,:,:,:,kc,ic)*t4r(:,:,:,:,kc,jc)

                  tmp3r(:,:,:,:,ic,jc) = tmp3r(:,:,:,:,ic,jc) +
     &               ux3r(:,:,:,:,ic,kc)*t3r(:,:,:,:,jc,kc) + ux3i(:,:,:,:,ic,kc)*t3i(:,:,:,:,jc,kc)

                  tmp3i(:,:,:,:,ic,jc) = tmp3i(:,:,:,:,ic,jc) -
     &               ux3r(:,:,:,:,ic,kc)*t3i(:,:,:,:,jc,kc) + ux3i(:,:,:,:,ic,kc)*t3r(:,:,:,:,jc,kc)

                  tmp2r(:,:,:,:,ic,jc) = tmp2r(:,:,:,:,ic,jc) +
     &               bx2r(:,:,:,:,ic,kc)*t2r(:,:,:,:,kc,jc) - bx2i(:,:,:,:,ic,kc)*t2i(:,:,:,:,kc,jc)

                  tmp2i(:,:,:,:,ic,jc) = tmp2i(:,:,:,:,ic,jc) +
     &               bx2r(:,:,:,:,ic,kc)*t2i(:,:,:,:,kc,jc) + bx2i(:,:,:,:,ic,kc)*t2r(:,:,:,:,kc,jc)

                  tmp1r(:,:,:,:,ic,jc) = tmp1r(:,:,:,:,ic,jc) +
     &               ux1r(:,:,:,:,kc,ic)*t1r(:,:,:,:,jc,kc) - ux1i(:,:,:,:,kc,ic)*t1i(:,:,:,:,jc,kc)

                  tmp1i(:,:,:,:,ic,jc) = tmp1i(:,:,:,:,ic,jc) -
     &               ux1r(:,:,:,:,kc,ic)*t1i(:,:,:,:,jc,kc) - ux1i(:,:,:,:,kc,ic)*t1r(:,:,:,:,jc,kc)

               end do
            end do
         end do

         s1r = 0.d0
         s1i = 0.d0
         s2r = 0.d0
         s2i = 0.d0
         s3r = 0.d0
         s3i = 0.d0
         s4r = 0.d0
         s4i = 0.d0

         do ic = 1, nc
            do jc = 1, nc
               do kc = 1, nc

                  s4r(:,:,:,:,ic,jc) = s4r(:,:,:,:,ic,jc) +
     &               tmp4r(:,:,:,:,ic,kc)*ux4r(:,:,:,:,kc,jc) - tmp4i(:,:,:,:,ic,kc)*ux4i(:,:,:,:,kc,jc)

                  s4i(:,:,:,:,ic,jc) = s4i(:,:,:,:,ic,jc) +
     &               tmp4r(:,:,:,:,ic,kc)*ux4i(:,:,:,:,kc,jc) + tmp4i(:,:,:,:,ic,kc)*ux4r(:,:,:,:,kc,jc)

                  s3r(:,:,:,:,ic,jc) = s3r(:,:,:,:,ic,jc) +
     &               tmp3r(:,:,:,:,ic,kc)*bx3r(:,:,:,:,jc,kc) + tmp3i(:,:,:,:,ic,kc)*bx3i(:,:,:,:,jc,kc)

                  s3i(:,:,:,:,ic,jc) = s3i(:,:,:,:,ic,jc) -
     &               tmp3r(:,:,:,:,ic,kc)*bx3i(:,:,:,:,jc,kc) + tmp3i(:,:,:,:,ic,kc)*bx3r(:,:,:,:,jc,kc)

                  s2r(:,:,:,:,ic,jc) = s2r(:,:,:,:,ic,jc) +
     &               tmp2r(:,:,:,:,ic,kc)*ux2r(:,:,:,:,jc,kc) + tmp2i(:,:,:,:,ic,kc)*ux2i(:,:,:,:,jc,kc)

                  s2i(:,:,:,:,ic,jc) = s2i(:,:,:,:,ic,jc) -
     &               tmp2r(:,:,:,:,ic,kc)*ux2i(:,:,:,:,jc,kc) + tmp2i(:,:,:,:,ic,kc)*ux2r(:,:,:,:,jc,kc)

                  s1r(:,:,:,:,ic,jc) = s1r(:,:,:,:,ic,jc) +
     &               tmp1r(:,:,:,:,ic,kc)*bx1r(:,:,:,:,kc,jc) - tmp1i(:,:,:,:,ic,kc)*bx1i(:,:,:,:,kc,jc)

                  s1i(:,:,:,:,ic,jc) = s1i(:,:,:,:,ic,jc) +
     &               tmp1r(:,:,:,:,ic,kc)*bx1i(:,:,:,:,kc,jc) + tmp1i(:,:,:,:,ic,kc)*bx1r(:,:,:,:,kc,jc)

               end do
            end do
         end do
c
c         ic(b),jc(u)
c        _|_ kc(b),mc(u)
c
         uar = cshift(bar, dim = that , shift =it)
         uai = cshift(bai, dim = that , shift =it)
         ubr = cshift(bbr, dim = that , shift =it)
         ubi = cshift(bbi, dim = that , shift =it)
         ucr = cshift(bcr, dim = that , shift =it)
         uci = cshift(bci, dim = that , shift =it)
         udr = cshift(bdr, dim = that , shift =it)
         udi = cshift(bdi, dim = that , shift =it)

         call link_conjug(bar,bai)
         call link_conjug(bbr,bbi)
         call link_conjug(ucr,uci)
         call link_conjug(udr,udi)

         do ix = 1, nx
            do iy = 1, ny
               do iz = 1, nz
                  do it = 1, nt

                     c1r = 0.d0
                     c1i = 0.d0

                     do ic = 1, nc
                       do jc = 1, nc
                         do kc = 1, nc
                           do mc = 1, nc

                              do ie = 1, 36

                                 c1r(ic,jc,kc,mc) = c1r(ic,jc,kc,mc) + fper(ie) * (
     &     ((( bcr(ix,iy,iz,it,lb(ie),kc)*bdr(ix,iy,iz,it,lc(ie),ic) -
     &         bci(ix,iy,iz,it,lb(ie),kc)*bdi(ix,iy,iz,it,lc(ie),ic))*s4r(ix,iy,iz,it,la(ie),lap(ie)) +
     &       (-bci(ix,iy,iz,it,lb(ie),kc)*bdr(ix,iy,iz,it,lc(ie),ic) -
     &         bcr(ix,iy,iz,it,lb(ie),kc)*bdi(ix,iy,iz,it,lc(ie),ic))*s4i(ix,iy,iz,it,la(ie),lap(ie)))*ucr(ix,iy,iz,it,mc,lbp(ie)) +
     &      ((-bci(ix,iy,iz,it,lb(ie),kc)*bdr(ix,iy,iz,it,lc(ie),ic) -
     &         bcr(ix,iy,iz,it,lb(ie),kc)*bdi(ix,iy,iz,it,lc(ie),ic))*s4r(ix,iy,iz,it,la(ie),lap(ie)) +
     &       ( bci(ix,iy,iz,it,lb(ie),kc)*bdi(ix,iy,iz,it,lc(ie),ic) -
     &         bcr(ix,iy,iz,it,lb(ie),kc)*bdr(ix,iy,iz,it,lc(ie),ic))*s4i(ix,iy,iz,it,la(ie),lap(ie)))*uci(ix,iy,iz,it,mc,lbp(ie)))*udr(ix,iy,iz,it,jc,lcp(ie)) +
     &     (((-bci(ix,iy,iz,it,lb(ie),kc)*bdr(ix,iy,iz,it,lc(ie),ic) -
     &         bcr(ix,iy,iz,it,lb(ie),kc)*bdi(ix,iy,iz,it,lc(ie),ic))*s4r(ix,iy,iz,it,la(ie),lap(ie)) +
     &       ( bci(ix,iy,iz,it,lb(ie),kc)*bdi(ix,iy,iz,it,lc(ie),ic) -
     &         bcr(ix,iy,iz,it,lb(ie),kc)*bdr(ix,iy,iz,it,lc(ie),ic))*s4i(ix,iy,iz,it,la(ie),lap(ie)))*ucr(ix,iy,iz,it,mc,lbp(ie)) +
     &      (( bci(ix,iy,iz,it,lb(ie),kc)*bdi(ix,iy,iz,it,lc(ie),ic) -
     &         bcr(ix,iy,iz,it,lb(ie),kc)*bdr(ix,iy,iz,it,lc(ie),ic))*s4r(ix,iy,iz,it,la(ie),lap(ie)) +
     &       ( bci(ix,iy,iz,it,lb(ie),kc)*bdr(ix,iy,iz,it,lc(ie),ic) +
     &         bcr(ix,iy,iz,it,lb(ie),kc)*bdi(ix,iy,iz,it,lc(ie),ic))*s4i(ix,iy,iz,it,la(ie),lap(ie)))*uci(ix,iy,iz,it,mc,lbp(ie)))*udi(ix,iy,iz,it,jc,lcp(ie))
     &                               )

                                  c1i(ic,jc,kc,mc) = c1i(ic,jc,kc,mc) + fper(ie) * (
     &     ((( bci(ix,iy,iz,it,lb(ie),kc)*bdr(ix,iy,iz,it,lc(ie),ic) +
     &         bcr(ix,iy,iz,it,lb(ie),kc)*bdi(ix,iy,iz,it,lc(ie),ic))*s4r(ix,iy,iz,it,la(ie),lap(ie)) +
     &       ( bcr(ix,iy,iz,it,lb(ie),kc)*bdr(ix,iy,iz,it,lc(ie),ic) -
     &         bci(ix,iy,iz,it,lb(ie),kc)*bdi(ix,iy,iz,it,lc(ie),ic))*s4i(ix,iy,iz,it,la(ie),lap(ie)))*ucr(ix,iy,iz,it,mc,lbp(ie)) +
     &      (( bcr(ix,iy,iz,it,lb(ie),kc)*bdr(ix,iy,iz,it,lc(ie),ic) -
     &         bci(ix,iy,iz,it,lb(ie),kc)*bdi(ix,iy,iz,it,lc(ie),ic))*s4r(ix,iy,iz,it,la(ie),lap(ie)) +
     &       (-bci(ix,iy,iz,it,lb(ie),kc)*bdr(ix,iy,iz,it,lc(ie),ic) -
     &         bcr(ix,iy,iz,it,lb(ie),kc)*bdi(ix,iy,iz,it,lc(ie),ic))*s4i(ix,iy,iz,it,la(ie),lap(ie)))*uci(ix,iy,iz,it,mc,lbp(ie)))*udr(ix,iy,iz,it,jc,lcp(ie)) +
     &     ((( bcr(ix,iy,iz,it,lb(ie),kc)*bdr(ix,iy,iz,it,lc(ie),ic) -
     &         bci(ix,iy,iz,it,lb(ie),kc)*bdi(ix,iy,iz,it,lc(ie),ic))*s4r(ix,iy,iz,it,la(ie),lap(ie)) +
     &       (-bci(ix,iy,iz,it,lb(ie),kc)*bdr(ix,iy,iz,it,lc(ie),ic) -
     &         bcr(ix,iy,iz,it,lb(ie),kc)*bdi(ix,iy,iz,it,lc(ie),ic))*s4i(ix,iy,iz,it,la(ie),lap(ie)))*ucr(ix,iy,iz,it,mc,lbp(ie)) +
     &      ((-bci(ix,iy,iz,it,lb(ie),kc)*bdr(ix,iy,iz,it,lc(ie),ic) -
     &         bcr(ix,iy,iz,it,lb(ie),kc)*bdi(ix,iy,iz,it,lc(ie),ic))*s4r(ix,iy,iz,it,la(ie),lap(ie)) +
     &       ( bci(ix,iy,iz,it,lb(ie),kc)*bdi(ix,iy,iz,it,lc(ie),ic) -
     &         bcr(ix,iy,iz,it,lb(ie),kc)*bdr(ix,iy,iz,it,lc(ie),ic))*s4i(ix,iy,iz,it,la(ie),lap(ie)))*uci(ix,iy,iz,it,mc,lbp(ie)))*udi(ix,iy,iz,it,jc,lcp(ie))
     &                               )

                              end do

                           end do
                         end do
                       end do
                     end do
c
c        _ _ ic(b),jc(u)
c        _|_
c            kc(b),mc(u)
c

                     c2r = 0.d0
                     c2i = 0.d0

                     do ic = 1, nc
                       do jc = 1, nc
                         do kc = 1, nc
                           do mc = 1, nc

                              do ie = 1, 36
                                 c2r(ic,jc,kc,mc) = c2r(ic,jc,kc,mc) + fper(ie) * (
     &      ((bar(ix,iy,iz,it,ic,lcp(ie))*c1r(lap(ie),la(ie),kc,mc) -
     &        bai(ix,iy,iz,it,ic,lcp(ie))*c1i(lap(ie),la(ie),kc,mc))*s1r(ix,iy,iz,it,lb(ie),lbp(ie)) +
     &      (-bai(ix,iy,iz,it,ic,lcp(ie))*c1r(lap(ie),la(ie),kc,mc) -
     &        bar(ix,iy,iz,it,ic,lcp(ie))*c1i(lap(ie),la(ie),kc,mc))*s1i(ix,iy,iz,it,lb(ie),lbp(ie)))*uar(ix,iy,iz,it,lc(ie),jc) +
     &     ((-bai(ix,iy,iz,it,ic,lcp(ie))*c1r(lap(ie),la(ie),kc,mc) -
     &        bar(ix,iy,iz,it,ic,lcp(ie))*c1i(lap(ie),la(ie),kc,mc))*s1r(ix,iy,iz,it,lb(ie),lbp(ie)) +
     &       (bai(ix,iy,iz,it,ic,lcp(ie))*c1i(lap(ie),la(ie),kc,mc) -
     &        bar(ix,iy,iz,it,ic,lcp(ie))*c1r(lap(ie),la(ie),kc,mc))*s1i(ix,iy,iz,it,lb(ie),lbp(ie)))*uai(ix,iy,iz,it,lc(ie),jc))

                                  c2i(ic,jc,kc,mc) = c2i(ic,jc,kc,mc) + fper(ie) * (
     &      ((bai(ix,iy,iz,it,ic,lcp(ie))*c1r(lap(ie),la(ie),kc,mc) +
     &        bar(ix,iy,iz,it,ic,lcp(ie))*c1i(lap(ie),la(ie),kc,mc))*s1r(ix,iy,iz,it,lb(ie),lbp(ie)) +
     &      ( bar(ix,iy,iz,it,ic,lcp(ie))*c1r(lap(ie),la(ie),kc,mc) -
     &        bai(ix,iy,iz,it,ic,lcp(ie))*c1i(lap(ie),la(ie),kc,mc))*s1i(ix,iy,iz,it,lb(ie),lbp(ie)))*uar(ix,iy,iz,it,lc(ie),jc) +
     &     (( bar(ix,iy,iz,it,ic,lcp(ie))*c1r(lap(ie),la(ie),kc,mc) -
     &        bai(ix,iy,iz,it,ic,lcp(ie))*c1i(lap(ie),la(ie),kc,mc))*s1r(ix,iy,iz,it,lb(ie),lbp(ie)) +
     &      (-bai(ix,iy,iz,it,ic,lcp(ie))*c1r(lap(ie),la(ie),kc,mc) -
     &        bar(ix,iy,iz,it,ic,lcp(ie))*c1i(lap(ie),la(ie),kc,mc))*s1i(ix,iy,iz,it,lb(ie),lbp(ie)))*uai(ix,iy,iz,it,lc(ie),jc))

                              end do

                           end do
                         end do
                       end do
                     end do

c
c        _ _ _
c        _|_| ic(b),jc(u)
c           kc(b),mc(u)
c

                     c1r = 0.d0
                     c1i = 0.d0

                     do ic = 1, nc
                       do jc = 1, nc
                         do kc = 1, nc
                           do mc = 1, nc

                              do ie = 1, 36
                                 c1r(ic,jc,kc,mc) = c1r(ic,jc,kc,mc) + fper(ie) * (
     &      ((bbr(ix,iy,iz,it,lb(ie),ic)*c2r(la(ie),lap(ie),kc,mc) -
     &        bbi(ix,iy,iz,it,lb(ie),ic)*c2i(la(ie),lap(ie),kc,mc))*s2r(ix,iy,iz,it,lb(ie),lbp(ie)) +
     &      (-bbi(ix,iy,iz,it,lb(ie),ic)*c2r(la(ie),lap(ie),kc,mc) -
     &        bbr(ix,iy,iz,it,lb(ie),ic)*c2i(la(ie),lap(ie),kc,mc))*s2i(ix,iy,iz,it,lb(ie),lbp(ie)))*ubr(ix,iy,iz,it,jc,lbp(ie)) +
     &     ((-bbi(ix,iy,iz,it,lb(ie),ic)*c2r(la(ie),lap(ie),kc,mc) -
     &        bbr(ix,iy,iz,it,lb(ie),ic)*c2i(la(ie),lap(ie),kc,mc))*s2r(ix,iy,iz,it,lb(ie),lbp(ie)) +
     &       (bbi(ix,iy,iz,it,lb(ie),ic)*c2i(la(ie),lap(ie),kc,mc) -
     &        bbr(ix,iy,iz,it,lb(ie),ic)*c2r(la(ie),lap(ie),kc,mc))*s2i(ix,iy,iz,it,lb(ie),lbp(ie)))*ubi(ix,iy,iz,it,jc,lbp(ie)))

                                  c1i(ic,jc,kc,mc) = c1i(ic,jc,kc,mc) + fper(ie) * (
     &      ((bbi(ix,iy,iz,it,lb(ie),ic)*c2r(la(ie),lap(ie),kc,mc) +
     &        bbr(ix,iy,iz,it,lb(ie),ic)*c2i(la(ie),lap(ie),kc,mc))*s2r(ix,iy,iz,it,lb(ie),lbp(ie)) +
     &      ( bbr(ix,iy,iz,it,lb(ie),ic)*c2r(la(ie),lap(ie),kc,mc) -
     &        bbi(ix,iy,iz,it,lb(ie),ic)*c2i(la(ie),lap(ie),kc,mc))*s2i(ix,iy,iz,it,lb(ie),lbp(ie)))*ubr(ix,iy,iz,it,jc,lbp(ie)) +
     &     (( bbr(ix,iy,iz,it,lb(ie),ic)*c2r(la(ie),lap(ie),kc,mc) -
     &        bbi(ix,iy,iz,it,lb(ie),ic)*c2i(la(ie),lap(ie),kc,mc))*s2r(ix,iy,iz,it,lb(ie),lbp(ie)) +
     &      (-bbi(ix,iy,iz,it,lb(ie),ic)*c2r(la(ie),lap(ie),kc,mc) -
     &        bbr(ix,iy,iz,it,lb(ie),ic)*c2i(la(ie),lap(ie),kc,mc))*s2i(ix,iy,iz,it,lb(ie),lbp(ie)))*ubi(ix,iy,iz,it,jc,lbp(ie)))

                              end do

                           end do
                         end do
                       end do
                     end do

                     res(ix,iy,iz,it) = 0.d0

                     do ie = 1, 36

                        res(ix,iy,iz,it) = res(ix,iy,iz,it) + fper(ie) * (
     &                     c1r(lap(ie),la(ie),lbp(ie),lb(ie))*s3r(ix,iy,iz,it,lc(ie),lcp(ie)) -
     &                     c1i(lap(ie),la(ie),lbp(ie),lb(ie))*s3i(ix,iy,iz,it,lc(ie),lcp(ie)) )

                     end do

                 end do
               end do
            end do
         end do

         W(:,:,:,:,1,lt) = res(:,:,:,:)

      end do ! end lt loop

      return

      end subroutine quad_loops

      END MODULE L_LOOPS
