c
c     ----------------------------------------------
c     CASE 7 - quarks at (11,0), (0,-6) and (0, 6)
c     ----------------------------------------------
c
c     base (bx)
c
      call product(ur,ui,ydir,cxlr,cxli,6)

      dxlr = cshift( cxlr, dim=ydir,shift=-6)
      dxli = cshift( cxli, dim=ydir,shift=-6)

      include 'basex.f'
c
c     top diagonal
c
      mx2r = mx1r
      mx2i = mx1i

      call product(ur,ui,xdir,mx1r,mx1i, 6)

      my2r = my1r
      my2i = my1i

      call baseup(mx1r,mx1i,mx2r,mx2i,my1r,my1i,my2r,my2i,6,3,6,11,xdir,ydir,bllr,blli)
      call link_conjug(bllr,blli)
c
c     bottom diagonal
c
      call basedo(mx1r,mx1i,mx2r,mx2i,my1r,my1i,my2r,my2i,-6,-3,6,11,xdir,ydir,brlr,brli)
      call link_conjug(brlr,brli)
c
c     time links shifting
c
      do it= 1, nL(that)

         t_3r(:,:,:,:,:,:) = cshift( tlr(:,:,:,:,:,:,it), dim=ydir,shift=-6)
         t_3i(:,:,:,:,:,:) = cshift( tli(:,:,:,:,:,:,it), dim=ydir,shift=-6)

         t_1r(:,:,:,:,:,:) = cshift( tlr(:,:,:,:,:,:,it), dim=ydir,shift= 6)
         t_1i(:,:,:,:,:,:) = cshift( tli(:,:,:,:,:,:,it), dim=ydir,shift= 6)

         t_2r(:,:,:,:,:,:) = cshift( tlr(:,:,:,:,:,:,it), dim=xdir,shift= 11)
         t_2i(:,:,:,:,:,:) = cshift( tli(:,:,:,:,:,:,it), dim=xdir,shift= 11)

         include 'res.f'

         W(:,:,:,:,7,it) = Res(:,:,:,:) / 6.0d0

      end do