c
c     ----------------------------------------------
c     CASE 2 - quarks at ( 3,0), (0,-2) and (0, 2)
c     ----------------------------------------------
c
c     lower diagonal (bl) [top upside down]
c
      bllr = cshift( cshift( bllr, dim=xdir,shift=-3), dim=ydir, shift=-2)
      blli = cshift( cshift( blli, dim=xdir,shift=-3), dim=ydir, shift=-2)
c
c     bottom diagonal (br)
c
      brlr = cshift( cshift( brlr, dim=xdir,shift=-3), dim=ydir, shift= 2)
      brli = cshift( cshift( brli, dim=xdir,shift=-3), dim=ydir, shift= 2)
c
c     time links shifting
c
      do it= 1, nL(that)

         t_dr(:,:,:,:,:,:) = cshift( tlr(:,:,:,:,:,:,it), dim=ydir,shift=-2)
         t_di(:,:,:,:,:,:) = cshift( tli(:,:,:,:,:,:,it), dim=ydir,shift=-2)

         t_ur(:,:,:,:,:,:) = cshift( tlr(:,:,:,:,:,:,it), dim=ydir,shift= 2)
         t_ui(:,:,:,:,:,:) = cshift( tli(:,:,:,:,:,:,it), dim=ydir,shift= 2)

         t_tr(:,:,:,:,:,:) = cshift( tlr(:,:,:,:,:,:,it), dim=xdir,shift=-3)
         t_ti(:,:,:,:,:,:) = cshift( tli(:,:,:,:,:,:,it), dim=xdir,shift=-3)

         include 'res_M.f'

         W_M(:,:,:,:,2,it) = Res(:,:,:,:) / 6.0d0

      end do
