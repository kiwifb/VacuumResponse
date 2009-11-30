c
c     ----------------------------------------------
c     CASE 4 - quarks at ( 5,0), (0,-3) and (0, 3)
c     ----------------------------------------------
c
c     top diagonal
c
      bllr = cshift( cshift( bllr, dim=xdir,shift=-5), dim=ydir, shift=-3)
      blli = cshift( cshift( blli, dim=xdir,shift=-5), dim=ydir, shift=-3)
c
c     bottom diagonal (br)
c
      brlr = cshift( cshift( brlr, dim=xdir,shift=-5), dim=ydir, shift= 3)
      brli = cshift( cshift( brli, dim=xdir,shift=-5), dim=ydir, shift= 3)
c
c     time links shifting
c
      do it= 1, nL(that)

         t_dr(:,:,:,:,:,:) = cshift( tlr(:,:,:,:,:,:,it), dim=ydir,shift=-3)
         t_di(:,:,:,:,:,:) = cshift( tli(:,:,:,:,:,:,it), dim=ydir,shift=-3)

         t_ur(:,:,:,:,:,:) = cshift( tlr(:,:,:,:,:,:,it), dim=ydir,shift= 3)
         t_ui(:,:,:,:,:,:) = cshift( tli(:,:,:,:,:,:,it), dim=ydir,shift= 3)

         t_tr(:,:,:,:,:,:) = cshift( tlr(:,:,:,:,:,:,it), dim=xdir,shift=-5)
         t_ti(:,:,:,:,:,:) = cshift( tli(:,:,:,:,:,:,it), dim=xdir,shift=-5)

         include 'res_M.f'

         W_M(:,:,:,:,4,it) = Res(:,:,:,:) / 6.0d0

      end do
