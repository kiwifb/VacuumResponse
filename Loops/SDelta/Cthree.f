c
c     ----------------------------------------------
c     CASE 3 - quarks at ( 4,0), (0,-2) and (0, 2)
c     ----------------------------------------------
c
c     base (bx)
c
!     Is the same as in the previous loop - no change.

      call product(ur,ui,xdir,cylr,cyli,4)
c
c     top diagonal
c
      dylr = cshift( cylr, dim=ydir,shift= 2)
      dyli = cshift( cyli, dim=ydir,shift= 2)

      tmp1r = cshift( cxlr, dim=xdir,shift= 4)
      tmp1i = cshift( cxli, dim=xdir,shift= 4)

      include 'baseup.f'
c
c     bottom diagonal
c
      dylr = cshift( cylr, dim=ydir,shift=-2)
      dyli = cshift( cyli, dim=ydir,shift=-2)

      tmp1r = cshift( dxlr, dim=xdir,shift= 4)
      tmp1i = cshift( dxli, dim=xdir,shift= 4)

      include 'basedo.f'
c
c     time links shifting
c
      do it= 1, nL(that)

         t_3r(:,:,:,:,:,:) = cshift( tlr(:,:,:,:,:,:,it), dim=ydir,shift=-2)
         t_3i(:,:,:,:,:,:) = cshift( tli(:,:,:,:,:,:,it), dim=ydir,shift=-2)

         t_1r(:,:,:,:,:,:) = cshift( tlr(:,:,:,:,:,:,it), dim=ydir,shift= 2)
         t_1i(:,:,:,:,:,:) = cshift( tli(:,:,:,:,:,:,it), dim=ydir,shift= 2)

         t_2r(:,:,:,:,:,:) = cshift( tlr(:,:,:,:,:,:,it), dim=xdir,shift= 4)
         t_2i(:,:,:,:,:,:) = cshift( tli(:,:,:,:,:,:,it), dim=xdir,shift= 4)

         include 'res.f'

         W(:,:,:,:,3,it) = Res(:,:,:,:) / 6.0d0

      end do
