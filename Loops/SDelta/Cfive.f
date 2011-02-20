c
c     ----------------------------------------------
c     CASE 5 - quarks at ( 7,0), (0,-4) and (0, 4)
c     ----------------------------------------------
c
c     base (bx)
c
      call product(ur,ui,ydir,cxlr,cxli,4)

      dxlr = cshift( cxlr, dim=ydir,shift=-4)
      dxli = cshift( cxli, dim=ydir,shift=-4)

      include 'basex.f'
c
c     top diagonal
c

      call product(ur,ui,xdir,mx1r,mx1i,1)
      call product(ur,ui,xdir,mx1r,mx1i,2)
      call product(ur,ui,xdir,mx1r,mx1i,3)
      mx2r = mx1r
      mx2i = mx1i
      call product(ur,ui,xdir,mx1r,mx1i,4)

      call product(ur,ui,ydir,my1r,my1i,1)
      call product(ur,ui,ydir,my1r,my1i,2)
      my2r = my1r
      my2i = my1i


      call baseup(mx1r,mx1i,mx2r,mx2i,my1r,my1i,my2r,my2i,4,2,4,7,xdir,ydir,bllr,blli)
      call link_conjug(bllr,blli)
c
c     bottom diagonal
c
      call basedo(mx1r,mx1i,mx2r,mx2i,my1r,my1i,my2r,my2i,-4,-2,4,7,xdir,ydir,brlr,brli)
      call link_conjug(brlr,brli)
c
c     time links shifting
c
      do it= 1, nL(that)

         t_3r(:,:,:,:,:,:) = cshift( tlr(:,:,:,:,:,:,it), dim=ydir,shift=-4)
         t_3i(:,:,:,:,:,:) = cshift( tli(:,:,:,:,:,:,it), dim=ydir,shift=-4)

         t_1r(:,:,:,:,:,:) = cshift( tlr(:,:,:,:,:,:,it), dim=ydir,shift= 4)
         t_1i(:,:,:,:,:,:) = cshift( tli(:,:,:,:,:,:,it), dim=ydir,shift= 4)

         t_2r(:,:,:,:,:,:) = cshift( tlr(:,:,:,:,:,:,it), dim=xdir,shift= 7)
         t_2i(:,:,:,:,:,:) = cshift( tli(:,:,:,:,:,:,it), dim=xdir,shift= 7)

         include 'res.f'

         W(:,:,:,:,5,it) = Res(:,:,:,:) / 6.0d0

      end do
