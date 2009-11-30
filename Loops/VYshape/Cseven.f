c
c     ----------------------------------------------
c     CASE 7 - quarks at (11,0), (0,-6) and (0, 6)
c     ----------------------------------------------
c
c     top diagonal
c
      bllr = up2x1r
      blli = up2x1i

      call multboxes(xdir,ydir,-2,1,bllr,blli,up2x1r,up2x1i)
      tmp1r = bllr
      tmp1i = blli

      call multboxes(xdir,ydir,-4,2,bllr,blli,up3x2r,up3x2i)

      call multboxes(xdir,ydir,-7,4,bllr,blli, tmp1r, tmp1i)

      bllr = cshift( bllr,dim=xdir, shift= 11)
      blli = cshift( blli,dim=xdir, shift= 11)
c
c     bottom diagonal
c
      brlr = do2x1r
      brli = do2x1i

      call multboxes(xdir,ydir,-2,-1,brlr,brli,do2x1r,do2x1i)
      tmp1r = brlr
      tmp1i = brli

      call multboxes(xdir,ydir,-4,-2,brlr,brli,do3x2r,do3x2i)

      call multboxes(xdir,ydir,-7,-4,brlr,brli, tmp1r, tmp1i)

      brlr = cshift( brlr,dim=xdir, shift= 11)
      brli = cshift( brli,dim=xdir, shift= 11)
c
c     time links shifting
c
      do it= 1, nL(that)

         t_dr(:,:,:,:,:,:) = cshift( tlr(:,:,:,:,:,:,it), dim=ydir,shift=-6)
         t_di(:,:,:,:,:,:) = cshift( tli(:,:,:,:,:,:,it), dim=ydir,shift=-6)

         t_ur(:,:,:,:,:,:) = cshift( tlr(:,:,:,:,:,:,it), dim=ydir,shift= 6)
         t_ui(:,:,:,:,:,:) = cshift( tli(:,:,:,:,:,:,it), dim=ydir,shift= 6)

         t_tr(:,:,:,:,:,:) = cshift( tlr(:,:,:,:,:,:,it), dim=xdir,shift= 11)
         t_ti(:,:,:,:,:,:) = cshift( tli(:,:,:,:,:,:,it), dim=xdir,shift= 11)

         include 'res.f'

         W(:,:,:,:,5,it) = Res(:,:,:,:) / 6.0d0

      end do