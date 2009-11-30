!
!  Parity Transform Routine to flip sign of Q leaving S invariant
!
!  Authors: Jianbo Zhang and Derek Leinweber
!
      MODULE GS_PARITYTRANSFORM

      CONTAINS

      subroutine ParityTransform(ur,ui)

      USE GS_LATTICESIZE

      implicit none

!     global variables

      integer,parameter                                       :: nc=3
      integer,parameter                                       :: mu=4
      double precision,dimension(nx,ny,nz,nt,mu,nc,nc)        :: ur,ui
!HPF$ DISTRIBUTE ur(*,*,BLOCK,BLOCK,*,*,*)
!HPF$ DISTRIBUTE ui(*,*,BLOCK,BLOCK,*,*,*)

!  Local Variables

      double precision,dimension(nx,ny,nz,nt,mu,nc,nc)        :: tr,ti
      integer :: iz, ipz, ipz1, iy, ipy, ipy1, ix, ipx, ipx1, ic, kc

!     parity transformation on U_mu(x)
!     P U_s(x) = U_s^\dagger(-x-\hat{s}),    PU_t(x) = U_t(-x)
!     U(nx,ny,nz,nt,mu,nc,nc)

      do iz = 1, nz
        ipz = -iz + nz
        if(ipz.eq.0) ipz =  nz
        ipz1 = ipz-1
        if(ipz1.eq.0) ipz1 =  nz

        do iy = 1, ny
          ipy = -iy + ny
          if(ipy.eq.0) ipy =  ny
          ipy1 = ipy-1
          if(ipy1.eq.0) ipy1 =  ny

          do ix = 1, nx
            ipx = -ix + nx
            if(ipx.eq.0)  ipx =  nx
            ipx1 = ipx-1
            if(ipx1.eq.0) ipx1 =  nx

            do kc = 1, nc
              do ic = 1, nc
                tr(ix,iy,iz,:,1,ic,kc) = ur(ipx1,ipy,ipz,:,1,kc,ic)
                ti(ix,iy,iz,:,1,ic,kc) =-ui(ipx1,ipy,ipz,:,1,kc,ic)

                tr(ix,iy,iz,:,2,ic,kc) = ur(ipx,ipy1,ipz,:,2,kc,ic)
                ti(ix,iy,iz,:,2,ic,kc) =-ui(ipx,ipy1,ipz,:,2,kc,ic)

                tr(ix,iy,iz,:,3,ic,kc) = ur(ipx,ipy,ipz1,:,3,kc,ic)
                ti(ix,iy,iz,:,3,ic,kc) =-ui(ipx,ipy,ipz1,:,3,kc,ic)
              end do
            end do

            tr(ix,iy,iz,:,4,:,:) = ur(ipx,ipy,ipz,:,4,:,:)
            ti(ix,iy,iz,:,4,:,:) = ui(ipx,ipy,ipz,:,4,:,:)

          end do
        end do
      end do

      ur = tr
      ui = ti

      return

      end subroutine parityTransform

      END MODULE GS_PARITYTRANSFORM
