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
      bllr = up2x1r
      blli = up2x1i

      call multboxes(xdir,ydir,-2,1,bllr,blli,up1x1r,up1x1i)

      call multboxes(xdir,ydir,-3,2,bllr,blli,up2x1r,up2x1i)

      bllr = cshift( bllr, dim=xdir, shift= 5)
      blli = cshift( blli, dim=xdir, shift= 5)
c
c     time links shifting
c
      do it= 1, nL(that)

         t_dr(:,:,:,:,:,:) = cshift( tlr(:,:,:,:,:,:,it), dim=ydir,shift=-3)
         t_di(:,:,:,:,:,:) = cshift( tli(:,:,:,:,:,:,it), dim=ydir,shift=-3)

         t_ur(:,:,:,:,:,:) = cshift( tlr(:,:,:,:,:,:,it), dim=ydir,shift= 3)
         t_ui(:,:,:,:,:,:) = cshift( tli(:,:,:,:,:,:,it), dim=ydir,shift= 3)

         t_tr(:,:,:,:,:,:) = cshift( tlr(:,:,:,:,:,:,it), dim=xdir,shift= 5)
         t_ti(:,:,:,:,:,:) = cshift( tli(:,:,:,:,:,:,it), dim=xdir,shift= 5)

         include 'res.f'

         W(:,:,:,:,4,it) = Res(:,:,:,:) / 6.0d0

      end do
