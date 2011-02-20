c
c     ----------------------------------------------
c     CASE 2 - quarks at ( 3,0), (0,-2) and (0, 2)
c     ----------------------------------------------
c
c     base (bx)
c
      call product(ur,ui,ydir,cxlr,cxli,3)
      call product(ur,ui,ydir,cxlr,cxli,4)

      bxlr = cshift( cxlr, dim=ydir,shift=-2)
      bxli = cshift( cxli, dim=ydir,shift=-2)
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

         t_3r(:,:,:,:,:,:) = cshift( tlr(:,:,:,:,:,:,it), dim=ydir,shift=-2)
         t_3i(:,:,:,:,:,:) = cshift( tli(:,:,:,:,:,:,it), dim=ydir,shift=-2)

         t_1r(:,:,:,:,:,:) = cshift( tlr(:,:,:,:,:,:,it), dim=ydir,shift= 2)
         t_1i(:,:,:,:,:,:) = cshift( tli(:,:,:,:,:,:,it), dim=ydir,shift= 2)

         t_2r(:,:,:,:,:,:) = cshift( tlr(:,:,:,:,:,:,it), dim=xdir,shift= 3)
         t_2i(:,:,:,:,:,:) = cshift( tli(:,:,:,:,:,:,it), dim=xdir,shift= 3)

         include 'res.f'

         W(:,:,:,:,2,it) = Res(:,:,:,:) / 6.0d0

      end do
