c
c     ----------------------------------------------
c     CASE 1 - quarks at ( 2,0), (0,-1) and (0, 1)
c     ----------------------------------------------
c
c     base (bx)
c
      call product(ur,ui,ydir,cxlr,cxli,1)
      call product(ur,ui,ydir,cxlr,cxli,2)

      bxlr = cshift( cxlr, dim=ydir,shift=-1)
      bxli = cshift( cxli, dim=ydir,shift=-1)
c
c     top diagonal (bl)
c
      bllr = cshift( up2x1r, dim=xdir,shift=2)
      blli = cshift( up2x1i, dim=xdir,shift=2)
c
c     time links shifting
c
      do it= 1, nL(that)

         t_ur(:,:,:,:,:,:) = cshift( tlr(:,:,:,:,:,:,it), dim=ydir,shift= 1)
         t_ui(:,:,:,:,:,:) = cshift( tli(:,:,:,:,:,:,it), dim=ydir,shift= 1)

         t_tr(:,:,:,:,:,:) = cshift( tlr(:,:,:,:,:,:,it), dim=xdir,shift= 2)
         t_ti(:,:,:,:,:,:) = cshift( tli(:,:,:,:,:,:,it), dim=xdir,shift= 2)

         t_dr(:,:,:,:,:,:) = cshift( tlr(:,:,:,:,:,:,it), dim=ydir,shift=-1)
         t_di(:,:,:,:,:,:) = cshift( tli(:,:,:,:,:,:,it), dim=ydir,shift=-1)

         include 'res.f'

         W(:,:,:,:,1,it) = Res(:,:,:,:) / 6.0d0

      end do
