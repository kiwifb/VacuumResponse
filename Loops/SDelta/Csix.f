c
c     ----------------------------------------------
c     CASE 6 - quarks at ( 9,0), (0,-5) and (0, 5)
c     ----------------------------------------------
c
c     base (bx)
c
      call product(ur,ui,ydir,cxlr,cxli, 5)

      dxlr = cshift( cxlr, dim=ydir,shift=-5)
      dxli = cshift( cxli, dim=ydir,shift=-5)

      include 'basex.f'
c
c     top diagonal
c
      mx2r = mx1r
      mx2i = mx1i

      call product(ur,ui,xdir,mx1r,mx1i, 5)

      call product(ur,ui,ydir,my1r,my1i, 3)

      call baseup(mx1r,mx1i,mx2r,mx2i,my1r,my1i,my2r,my2i,5,2,5,9,xdir,ydir,bllr,blli)
      call link_conjug(bllr,blli)
c
c     bottom diagonal
c
      call basedo(mx1r,mx1i,mx2r,mx2i,my1r,my1i,my2r,my2i,-5,-2,5,9,xdir,ydir,brlr,brli)
      call link_conjug(brlr,brli)
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
