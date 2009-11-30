c
c     ----------------------------------------------
c     CASE 2 - quarks at ( 3,0), (0,-2) and (0, 2)
c     ----------------------------------------------
c
c     top diagonal
c
      bllr = cshift( up3x2r, dim=xdir,shift= 3)
      blli = cshift( up3x2i, dim=xdir,shift= 3)
c
c     bottom diagonal
c
      brlr = cshift( do3x2r, dim=xdir,shift= 3)
      brli = cshift( do3x2i, dim=xdir,shift= 3)
c
c     time links shifting
c
      do it= 1, nL(that)

         t_dr(:,:,:,:,:,:) = cshift( tlr(:,:,:,:,:,:,it), dim=ydir,shift=-2)
         t_di(:,:,:,:,:,:) = cshift( tli(:,:,:,:,:,:,it), dim=ydir,shift=-2)

         t_ur(:,:,:,:,:,:) = cshift( tlr(:,:,:,:,:,:,it), dim=ydir,shift= 2)
         t_ui(:,:,:,:,:,:) = cshift( tli(:,:,:,:,:,:,it), dim=ydir,shift= 2)

         t_tr(:,:,:,:,:,:) = cshift( tlr(:,:,:,:,:,:,it), dim=xdir,shift= 3)
         t_ti(:,:,:,:,:,:) = cshift( tli(:,:,:,:,:,:,it), dim=xdir,shift= 3)

         include 'res.f'

         W(:,:,:,:,2,it) = Res(:,:,:,:) / 6.0d0

      end do
