c
c     ----------------------------------------------
c     CASE 6 - quarks at ( 9,0), (0,-5) and (0, 5)
c     ----------------------------------------------
c
c     base (bx)
c
      call product(ur,ui,ydir,cxlr,cxli, 9)
      call product(ur,ui,ydir,cxlr,cxli,10)

      bxlr = cshift( cxlr, dim=ydir,shift=-5)
      bxli = cshift( cxli, dim=ydir,shift=-5)
c
c     top diagonal
c
      call multboxes(xdir,ydir,-6, 4,ilr,ili,up2x1r,up2x1i)

      bllr = ilr
      blli = ili

      call multboxes(xdir,ydir,-8, 5,bllr,blli,ixr,ixi)

      bllr = cshift( bllr, dim=xdir, shift= 9)
      blli = cshift( blli, dim=xdir, shift= 9)
c
c     bottom diagonal
c
      call multboxes(xdir,ydir,-6,-4,irr,iri,do2x1r,do2x1i)

      brlr = irr
      brli = iri

      call multboxes(xdir,ydir,-8,-5,brlr,brli,ixr,ixi)

      brlr = cshift( brlr, dim=xdir, shift= 9)
      brli = cshift( brli, dim=xdir, shift= 9)
c
c     time links shifting
c
      do it= 1, nL(that)

         t_3r(:,:,:,:,:,:) = cshift( tlr(:,:,:,:,:,:,it), dim=ydir,shift=-5)
         t_3i(:,:,:,:,:,:) = cshift( tli(:,:,:,:,:,:,it), dim=ydir,shift=-5)

         t_1r(:,:,:,:,:,:) = cshift( tlr(:,:,:,:,:,:,it), dim=ydir,shift= 5)
         t_1i(:,:,:,:,:,:) = cshift( tli(:,:,:,:,:,:,it), dim=ydir,shift= 5)

         t_2r(:,:,:,:,:,:) = cshift( tlr(:,:,:,:,:,:,it), dim=xdir,shift= 9)
         t_2i(:,:,:,:,:,:) = cshift( tli(:,:,:,:,:,:,it), dim=xdir,shift= 9)

         include 'res.f'

         W(:,:,:,:,6,it) = Res(:,:,:,:) / 6.0d0

      end do
