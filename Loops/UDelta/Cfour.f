c
c     ----------------------------------------------
c     CASE 4 - quarks at ( 5,0), (0,-3) and (0, 3)
c     ----------------------------------------------
c
c     base (bx)
c
      call product(ur,ui,ydir,cxlr,cxli,3)

      dxlr = cshift( cxlr, dim=ydir,shift=-3)
      dxli = cshift( cxli, dim=ydir,shift=-3)

      include 'basex.f'

      call product(ur,ui,xdir,cylr,cyli,5)
c
c     top diagonal
c
      dylr = cshift( cylr, dim=ydir,shift= 3)
      dyli = cshift( cyli, dim=ydir,shift= 3)

      tmp1r = cshift( cxlr, dim=xdir,shift= 5)
      tmp1i = cshift( cxli, dim=xdir,shift= 5)

      include 'baseup.f'
c
c     bottom diagonal
c
      dylr = cshift( cylr, dim=ydir,shift=-3)
      dyli = cshift( cyli, dim=ydir,shift=-3)

      tmp1r = cshift( dxlr, dim=xdir,shift= 5)
      tmp1i = cshift( dxli, dim=xdir,shift= 5)

      include 'basedo.f'
c
c     time links shifting
c
      do it= 1, nL(that)

         t_3r(:,:,:,:,:,:) = cshift( tlr(:,:,:,:,:,:,it), dim=ydir,shift=-3)
         t_3i(:,:,:,:,:,:) = cshift( tli(:,:,:,:,:,:,it), dim=ydir,shift=-3)

         t_1r(:,:,:,:,:,:) = cshift( tlr(:,:,:,:,:,:,it), dim=ydir,shift= 3)
         t_1i(:,:,:,:,:,:) = cshift( tli(:,:,:,:,:,:,it), dim=ydir,shift= 3)

         t_2r(:,:,:,:,:,:) = cshift( tlr(:,:,:,:,:,:,it), dim=xdir,shift= 5)
         t_2i(:,:,:,:,:,:) = cshift( tli(:,:,:,:,:,:,it), dim=xdir,shift= 5)

         include 'res.f'

         W(:,:,:,:,4,it) = Res(:,:,:,:) / 6.0d0

      end do
