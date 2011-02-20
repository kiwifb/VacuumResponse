c
c     ----------------------------------------------
c     CASE 5 - quarks at ( 7,0), (0,-4) and (0, 4)
c     ----------------------------------------------
c
c     base (bx)
c
      call product(ur,ui,ydir,cxlr,cxli,7)
      call product(ur,ui,ydir,cxlr,cxli,8)

      bxlr = cshift( cxlr, dim=ydir,shift=-4)
      bxli = cshift( cxli, dim=ydir,shift=-4)
c
c     top diagonal
c
      call multboxes(xdir,ydir,-4, 3,ilr,ili,up2x1r,up2x1i)

      bllr = ilr
      blli = ili

      call multboxes(xdir,ydir,-6, 4,bllr,blli,ixr,ixi)

      bllr = cshift( bllr, dim=xdir, shift= 7)
      blli = cshift( blli, dim=xdir, shift= 7)
c
c     bottom diagonal
c
      call multboxes(xdir,ydir,-4,-3,irr,iri,do2x1r,do2x1i)

      brlr = irr
      brli = iri

      call multboxes(xdir,ydir,-6,-4,brlr,brli,ixr,ixi)

      brlr = cshift( brlr, dim=xdir, shift= 7)
      brli = cshift( brli, dim=xdir, shift= 7)
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
