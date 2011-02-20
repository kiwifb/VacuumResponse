c
c     ----------------------------------------------
c     CASE 4 - quarks at ( 5,0), (0,-3) and (0, 3)
c     ----------------------------------------------
c
c     base (bx)
c
      call product(ur,ui,ydir,cxlr,cxli,5)
      call product(ur,ui,ydir,cxlr,cxli,6)

      bxlr = cshift( cxlr, dim=ydir,shift=-3)
      bxli = cshift( cxli, dim=ydir,shift=-3)
c
c     top diagonal
c
      ilr = uyr
      ili = uyi

      call multboxes(xdir,ydir, 0, 1,ilr,ili,up2x1r,up2x1i)
      call multboxes(xdir,ydir,-2, 2,ilr,ili,up2x1r,up2x1i)

      bllr = ilr
      blli = ili

      call multboxes(xdir,ydir,-4, 3,bllr,blli,ixr,ixi)

      bllr = cshift( bllr, dim=xdir, shift= 5)
      blli = cshift( blli, dim=xdir, shift= 5)
c
c     bottom diagonal
c
      irr = dyr
      iri = dyi

      call multboxes(xdir,ydir, 0,-1,irr,iri,do2x1r,do2x1i)
      call multboxes(xdir,ydir,-2,-2,irr,iri,do2x1r,do2x1i)

      brlr = irr
      brli = iri

      call multboxes(xdir,ydir,-4,-3,brlr,brli,ixr,ixi)

      brlr = cshift( brlr, dim=xdir, shift= 5)
      brli = cshift( brli, dim=xdir, shift= 5)
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
