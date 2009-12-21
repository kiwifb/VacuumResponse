c
c     ----------------------------------------------
c     CASE 7 - quarks at (11,0), (0,-6) and (0, 6)
c     ----------------------------------------------
c
c     base (bx)
c
      call product(ur,ui,ydir,cxlr,cxli,6)

      dxlr = cshift( cxlr, dim=ydir,shift=-6)
      dxli = cshift( cxli, dim=ydir,shift=-6)

      include 'basex.f'

      call product(ur,ui,xdir,cylr,cyli,10)
      call product(ur,ui,xdir,cylr,cyli,11)
c
c     top diagonal
c
      dylr = cshift( cylr, dim=ydir,shift= 6)
      dyli = cshift( cyli, dim=ydir,shift= 6)

      tmp1r = cshift( cxlr, dim=xdir,shift=11)
      tmp1i = cshift( cxli, dim=xdir,shift=11)

      include 'baseup.f'
c
c     bottom diagonal
c
      dylr = cshift( cylr, dim=ydir,shift=-6)
      dyli = cshift( cyli, dim=ydir,shift=-6)

      tmp1r = cshift( dxlr, dim=xdir,shift=11)
      tmp1i = cshift( dxli, dim=xdir,shift=11)

      include 'basedo.f'
c
c     time links shifting
c
      do it= 1, nL(that)

         t_3r(:,:,:,:,:,:) = cshift( tlr(:,:,:,:,:,:,it), dim=ydir,shift=-6)
         t_3i(:,:,:,:,:,:) = cshift( tli(:,:,:,:,:,:,it), dim=ydir,shift=-6)

         t_1r(:,:,:,:,:,:) = cshift( tlr(:,:,:,:,:,:,it), dim=ydir,shift= 6)
         t_1i(:,:,:,:,:,:) = cshift( tli(:,:,:,:,:,:,it), dim=ydir,shift= 6)

         t_2r(:,:,:,:,:,:) = cshift( tlr(:,:,:,:,:,:,it), dim=xdir,shift= 11)
         t_2i(:,:,:,:,:,:) = cshift( tli(:,:,:,:,:,:,it), dim=xdir,shift= 11)

         include 'res.f'

         W(:,:,:,:,7,it) = Res(:,:,:,:) / 6.0d0

      end do