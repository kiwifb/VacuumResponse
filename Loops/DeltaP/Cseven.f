c
c     ----------------------------------------------
c     CASE 7 - quarks at (11,0), (0,-6) and (0, 6)
c     ----------------------------------------------
c
c     base (bx)
c
      call product(ur,ui,ydir,cxlr,cxli,11)
      call product(ur,ui,ydir,cxlr,cxli,12)

      bxlr = cshift( cxlr, dim=ydir,shift=-6)
      bxli = cshift( cxli, dim=ydir,shift=-6)
c
c     top diagonal
c
      call multboxes(xdir,ydir,-8, 5,ilr,ili,up2x1r,up2x1i)

      bllr = ilr
      blli = ili

      call multboxes(xdir,ydir,-10, 6,bllr,blli,ixr,ixi)

      bllr = cshift( bllr,dim=xdir, shift= 11)
      blli = cshift( blli,dim=xdir, shift= 11)
c
c     bottom diagonal
c
      call multboxes(xdir,ydir,-8,-5,irr,iri,do2x1r,do2x1i)

      brlr = irr
      brli = iri

      call multboxes(xdir,ydir,-10,-6,brlr,brli,ixr,ixi)

      brlr = cshift( brlr,dim=xdir, shift= 11)
      brli = cshift( brli,dim=xdir, shift= 11)
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