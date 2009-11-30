c
c     ----------------------------------------------
c     CASE 2 - quarks at ( 3,0), (0,-2) and (0, 2)
c     ----------------------------------------------
c
c     base
c
      call product(ur,ui,ydir,bxlr,bxli,3)
      call product(ur,ui,ydir,bxlr,bxli,4)
c
c     top diagonal
c
      bllr = cshift( cshift( ld3r, dim=ydir,shift= 2),dim=xdir,shift= 3)
      blli = cshift( cshift( ld3i, dim=ydir,shift= 2),dim=xdir,shift= 3)
c
c     bottom diagonal
c
      brlr = cshift( cshift( rd3r, dim=ydir,shift= 2),dim=xdir,shift= 3)
      brli = cshift( cshift( rd3i, dim=ydir,shift= 2),dim=xdir,shift= 3)
c
c     time links shifting
c
      do it = 1, nL(that)

         tlr(:,:,:,:,:,:) = tltr(:,:,:,:,:,:,it)
         tli(:,:,:,:,:,:) = tlti(:,:,:,:,:,:,it)

         tl1r = cshift( tlr, dim=ydir,shift= 4)
         tl1i = cshift( tli, dim=ydir,shift= 4)

         tl2r = cshift( cshift( tlr, dim=ydir,shift= 2),dim=xdir,shift= 3)
         tl2i = cshift( cshift( tli, dim=ydir,shift= 2),dim=xdir,shift= 3)

         include 'YLfiles/res.f'

         Wp(:,:,:,:,2,it) = Res(:,:,:,:) / 216.0d0

         call delta_mirror(xdir,ydir,bxlr,bxli,bllr,blli,brlr,brli,tlr,tli,tl1r,tl1i,Res,it,3,2)

         Wm(:,:,:,:,2,it) = Res(:,:,:,:) / 216.0d0

      end do