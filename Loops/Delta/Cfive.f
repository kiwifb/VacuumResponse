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
      bllr = up2x1r
      blli = up2x1i

      call multboxes(xdir,ydir,-2,1,bllr,blli,up3x2r,up3x2i)

      call multboxes(xdir,ydir,-5,3,bllr,blli,up2x1r,up2x1i)

      bllr = cshift( bllr, dim=xdir, shift= 7)
      blli = cshift( blli, dim=xdir, shift= 7)
c
c     bottom diagonal
c
      brlr = do2x1r
      brli = do2x1i

      call multboxes(xdir,ydir,-2,-1,brlr,brli,do3x2r,do3x2i)

      call multboxes(xdir,ydir,-5,-3,brlr,brli,do2x1r,do2x1i)

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
