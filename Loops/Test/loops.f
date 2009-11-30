c
c     subroutine to calculate valules for Wilson loops, W, for 'VY' shape
!        adapted from 'Delta' shape
c     1.0 20060110 FB
!     2.0 FB071129 - ported to use boxes.f and hopefully correct a bug causing gauge
!        invariance problem in shape 7.
c
      MODULE L_LOOPS

      CONTAINS

c========================================================================================================

      subroutine test_loops(xdir,ydir,ur,ui,W)
c
c     modules
c
      USE L_baryonParam
      USE L_product
      USE L_BOXES
c
c     no implicit typing
c
      implicit none
c
c     global variables
c
#include"loopsize.f"
      integer                                                           :: xdir,ydir
      double precision,dimension(nx,ny,nz,nt,mu,nc,nc)            	 	:: ur,ui
!HPF$ DISTRIBUTE ur(*,*,BLOCK,BLOCK,*,*,*)
!HPF$ DISTRIBUTE ui(*,*,BLOCK,BLOCK,*,*,*)
c
c     the five cases (1-nYloop) are stored in the second to last place in 'W'
c
      double precision,dimension(nx,ny,nz,nt,0:nloop,nL(that))          :: W
!HPF$ DISTRIBUTE W(*,*,BLOCK,BLOCK,*,*)
c
      include'loop_declarations.f'
c
c     ----------------
c     execution begins
c     ----------------
c
c     -----------------------
c     make the diagonal links
c     -----------------------
c
!      call boxm1p1(xdir,ydir,ur,ui,up1x1r,up1x1i)
      call boxm1m1(xdir,ydir,ur,ui,do1x1r,do1x1i)
!      call boxm2p1(xdir,ydir,ur,ui,up2x1r,up2x1i)
!      call boxm2m1(xdir,ydir,ur,ui,do2x1r,do2x1i)
!      call boxm3p2(xdir,ydir,ur,ui,up3x2r,up3x2i)
!      call boxm3m2(xdir,ydir,ur,ui,do3x2r,do3x2i)
c
c     ---------------------------------
c     start of Wilson loop calculations
c     ---------------------------------
c
c     compute the time link loop
c
      do it = 1, nL(that)

         call product(ur,ui,that,tmp1r,tmp1i,it)
         tlr(:,:,:,:,:,:,it) = tmp1r(:,:,:,:,:,:)
         tli(:,:,:,:,:,:,it) = tmp1i(:,:,:,:,:,:)

      end do

      brlr = cshift( do1x1r, dim=xdir,shift=1)
      brli = cshift( do1x1i, dim=xdir,shift=1)

      do it= 1, nL(that)

         t_ur(:,:,:,:,:,:) = cshift( tlr(:,:,:,:,:,:,it), dim=ydir,shift=-1)
         t_ui(:,:,:,:,:,:) = cshift( tli(:,:,:,:,:,:,it), dim=ydir,shift=-1)

         t_tr(:,:,:,:,:,:) = cshift( tlr(:,:,:,:,:,:,it), dim=xdir,shift= 1)
         t_ti(:,:,:,:,:,:) = cshift( tli(:,:,:,:,:,:,it), dim=xdir,shift= 1)

         sbrlr = cshift(brlr,dim=that,shift=it)
         sbrli = cshift(brli,dim=that,shift=it)

         tmp1r = 0.0d0
         tmp1i = 0.0d0

         tmp2r = 0.0d0
         tmp2i = 0.0d0

         do ic = 1,nc
            do jc = 1,nc
               do kc = 1,nc

                  tmp1r(:,:,:,:,ic,jc) = tmp1r(:,:,:,:,ic,jc)
     &                 + brlr(:,:,:,:,ic,kc) * t_ur(:,:,:,:,kc,jc) - brli(:,:,:,:,ic,kc) * t_ui(:,:,:,:,kc,jc)

                  tmp1i(:,:,:,:,ic,jc) = tmp1i(:,:,:,:,ic,jc)
     &                 + brlr(:,:,:,:,ic,kc) * t_ui(:,:,:,:,kc,jc) + brli(:,:,:,:,ic,kc) * t_ur(:,:,:,:,kc,jc)

                  tmp2r(:,:,:,:,ic,jc) = tmp2r(:,:,:,:,ic,jc)
     &                 + t_tr(:,:,:,:,ic,kc) * sbrlr(:,:,:,:,kc,jc) - t_ti(:,:,:,:,ic,kc) * sbrli(:,:,:,:,kc,jc)

                  tmp2i(:,:,:,:,ic,jc) = tmp2i(:,:,:,:,ic,jc)
     &                 + t_tr(:,:,:,:,ic,kc) * sbrli(:,:,:,:,kc,jc) + t_ti(:,:,:,:,ic,kc) * sbrlr(:,:,:,:,kc,jc)

               end do
            end do
         end do

         res = 0.d0

         do ic = 1, nc
            do jc = 1, nc

               res(:,:,:,:) = res(:,:,:,:) +
     &            tmp1r(:,:,:,:,ic,jc)*tmp2r(:,:,:,:,ic,jc) + tmp1i(:,:,:,:,ic,jc)*tmp2i(:,:,:,:,ic,jc)

            end do
         end do

         W(:,:,:,:,1,it) = Res(:,:,:,:)

      end do


      return

      end subroutine test_loops

      END MODULE L_LOOPS
