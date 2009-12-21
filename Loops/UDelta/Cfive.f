c
c     ----------------------------------------------
c     CASE 5 - quarks at ( 7,0), (0,-4) and (0, 4)
c     ----------------------------------------------
c
c     base (bx)
c
      call product(ur,ui,ydir,cxlr,cxli,4)

      dxlr = cshift( cxlr, dim=ydir,shift=-4)
      dxli = cshift( cxli, dim=ydir,shift=-4)

      include 'basex.f'

      call product(ur,ui,xdir,cylr,cyli,6)
      call product(ur,ui,xdir,cylr,cyli,7)
c
c     top diagonal
c
      dylr = cshift( cylr, dim=ydir,shift= 4)
      dyli = cshift( cyli, dim=ydir,shift= 4)

      tmp1r = cshift( cxlr, dim=xdir,shift= 7)
      tmp1i = cshift( cxli, dim=xdir,shift= 7)

      include 'baseup.f'
c
c     bottom diagonal
c
      dylr = cshift( cylr, dim=ydir,shift=-4)
      dyli = cshift( cyli, dim=ydir,shift=-4)

      tmp1r = cshift( dxlr, dim=xdir,shift= 7)
      tmp1i = cshift( dxli, dim=xdir,shift= 7)

      include 'basedo.f'
c
c     time links shifting
c
      do it= 1, nL(that)

         t_3r(:,:,:,:,:,:) = cshift( tlr(:,:,:,:,:,:,it), dim=ydir,shift=-4)
         t_3i(:,:,:,:,:,:) = cshift( tli(:,:,:,:,:,:,it), dim=ydir,shift=-4)

         t_1r(:,:,:,:,:,:) = cshift( tlr(:,:,:,:,:,:,it), dim=ydir,shift= 4)
         t_1i(:,:,:,:,:,:) = cshift( tli(:,:,:,:,:,:,it), dim=ydir,shift= 4)

         t_2r(:,:,:,:,:,:) = cshift( tlr(:,:,:,:,:,:,it), dim=xdir,shift= 7)
         t_2i(:,:,:,:,:,:) = cshift( tli(:,:,:,:,:,:,it), dim=xdir,shift= 7)

         include 'res.f'

         W(:,:,:,:,5,it) = Res(:,:,:,:) / 6.0d0

      end do
