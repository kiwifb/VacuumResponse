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
      bllr = up2x1r
      blli = up2x1i

      call multboxes(xdir,ydir,-2,1,bllr,blli,up2x1r,up2x1i)

      call multboxes(xdir,ydir,-4,2,bllr,blli,up1x1r,up1x1i)

      call multboxes(xdir,ydir,-5,3,bllr,blli,up2x1r,up2x1i)

      call multboxes(xdir,ydir,-7,4,bllr,blli,up2x1r,up2x1i)

      bllr = cshift( bllr, dim=xdir, shift= 9)
      blli = cshift( blli, dim=xdir, shift= 9)
c
c     bottom diagonal
c
      brlr = do2x1r
      brli = do2x1i

      call multboxes(xdir,ydir,-2,-1,brlr,brli,do2x1r,do2x1i)

      call multboxes(xdir,ydir,-4,-2,brlr,brli,do1x1r,do1x1i)

      call multboxes(xdir,ydir,-5,-3,brlr,brli,do2x1r,do2x1i)

      call multboxes(xdir,ydir,-7,-4,brlr,brli,do2x1r,do2x1i)

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
