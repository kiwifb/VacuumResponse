c
c     ----------------------------------------------
c     Mirror CASE 1 - quarks at ( -2,0), (0,-1) and (0, 1)
c     ----------------------------------------------
c
c     lower diagonal (bl) [top upside down]
c
      bllr = cshift( cshift( bllr, dim=xdir,shift=-2), dim=ydir, shift=-1)
      blli = cshift( cshift( blli, dim=xdir,shift=-2), dim=ydir, shift=-1)
c
c     bottom diagonal (br)
c
      brlr = cshift( cshift( brlr, dim=xdir,shift=-2), dim=ydir, shift= 1)
      brli = cshift( cshift( brli, dim=xdir,shift=-2), dim=ydir, shift= 1)
c
c     time links shifting
c
      do it= 1, nL(that)

         t_ur(:,:,:,:,:,:) = cshift( tlr(:,:,:,:,:,:,it), dim=ydir,shift= 1)
         t_ui(:,:,:,:,:,:) = cshift( tli(:,:,:,:,:,:,it), dim=ydir,shift= 1)

         t_tr(:,:,:,:,:,:) = cshift( tlr(:,:,:,:,:,:,it), dim=xdir,shift=-2)
         t_ti(:,:,:,:,:,:) = cshift( tli(:,:,:,:,:,:,it), dim=xdir,shift=-2)

         t_dr(:,:,:,:,:,:) = cshift( tlr(:,:,:,:,:,:,it), dim=ydir,shift=-1)
         t_di(:,:,:,:,:,:) = cshift( tli(:,:,:,:,:,:,it), dim=ydir,shift=-1)

         include 'res_M.f'

         W_M(:,:,:,:,1,it) = Res(:,:,:,:) / 6.0d0

      end do
