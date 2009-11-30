c
c     ----------------------------------------------
c     CASE 3 - quarks at ( 4,0), (0,-2) and (0, 2)
c     ----------------------------------------------
c
c     top diagonal
c
      bllr = up2x1r
      blli = up2x1i

      call multboxes(xdir,ydir,-2,1,bllr,blli,up2x1r,up2x1i)

      bllr = cshift( bllr, dim=xdir, shift= 4)
      blli = cshift( blli, dim=xdir, shift= 4)
c
c     bottom diagonal
c
      brlr = do2x1r
      brli = do2x1i

      call multboxes(xdir,ydir,-2,-1,brlr,brli,do2x1r,do2x1i)

      brlr = cshift( brlr, dim=xdir, shift= 4)
      brli = cshift( brli, dim=xdir, shift= 4)
c
c     time links shifting
c
      do it= 1, nL(that)

         t_dr(:,:,:,:,:,:) = cshift( tlr(:,:,:,:,:,:,it), dim=ydir,shift=-2)
         t_di(:,:,:,:,:,:) = cshift( tli(:,:,:,:,:,:,it), dim=ydir,shift=-2)

         t_ur(:,:,:,:,:,:) = cshift( tlr(:,:,:,:,:,:,it), dim=ydir,shift= 2)
         t_ui(:,:,:,:,:,:) = cshift( tli(:,:,:,:,:,:,it), dim=ydir,shift= 2)

         t_tr(:,:,:,:,:,:) = cshift( tlr(:,:,:,:,:,:,it), dim=xdir,shift= 4)
         t_ti(:,:,:,:,:,:) = cshift( tli(:,:,:,:,:,:,it), dim=xdir,shift= 4)

         include 'res.f'

         W(:,:,:,:,3,it) = Res(:,:,:,:) / 6.0d0

      end do
