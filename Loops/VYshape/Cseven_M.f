c
c     ----------------------------------------------
c     CASE 7 - quarks at (11,0), (0,-6) and (0, 6)
c     ----------------------------------------------
c
c     top diagonal
c
      bllr = cshift( cshift( bllr, dim=xdir,shift=-11), dim=ydir, shift=-6)
      blli = cshift( cshift( blli, dim=xdir,shift=-11), dim=ydir, shift=-6)
c
c     bottom diagonal (br)
c
      brlr = cshift( cshift( brlr, dim=xdir,shift=-11), dim=ydir, shift= 6)
      brli = cshift( cshift( brli, dim=xdir,shift=-11), dim=ydir, shift= 6)
c
c     time links shifting
c
      do it= 1, nL(that)

         t_dr(:,:,:,:,:,:) = cshift( tlr(:,:,:,:,:,:,it), dim=ydir,shift=-6)
         t_di(:,:,:,:,:,:) = cshift( tli(:,:,:,:,:,:,it), dim=ydir,shift=-6)

         t_ur(:,:,:,:,:,:) = cshift( tlr(:,:,:,:,:,:,it), dim=ydir,shift= 6)
         t_ui(:,:,:,:,:,:) = cshift( tli(:,:,:,:,:,:,it), dim=ydir,shift= 6)

         t_tr(:,:,:,:,:,:) = cshift( tlr(:,:,:,:,:,:,it), dim=xdir,shift=-11)
         t_ti(:,:,:,:,:,:) = cshift( tli(:,:,:,:,:,:,it), dim=xdir,shift=-11)

         include 'res_M.f'

         W_M(:,:,:,:,5,it) = Res(:,:,:,:) / 6.0d0

      end do