c
c     ----------------------------------------------
c     CASE 1 - quarks at ( 2,0), (0,-1) and (0, 1)
c     ----------------------------------------------
c
c     top diagonal (bl)
c
      bllr = cshift( ld2r,dim=xdir,shift=2)
      blli = cshift( ld2i,dim=xdir,shift=2)
c
c     bottom diagonal (br)
c
      brlr = cshift( rd2r,dim=xdir,shift=2)
      brli = cshift( rd2i,dim=xdir,shift=2)
c
c     time links shifting
c
      do it = 1, nL(that)

         tlr(:,:,:,:,:,:) = tltr(:,:,:,:,:,:,it)
         tli(:,:,:,:,:,:) = tlti(:,:,:,:,:,:,it)

         tl1r = cshift( tlr, dim=ydir,shift= 1)
         tl1i = cshift( tli, dim=ydir,shift= 1)

         tl2r = cshift( tlr, dim=xdir,shift= 2)
         tl2i = cshift( tli, dim=xdir,shift= 2)

         tlr = cshift( tlr, dim=ydir,shift=-1)
         tli = cshift( tli, dim=ydir,shift=-1)

         include 'YLfiles/res.f'

         W(:,:,:,:,1,it) = Res(:,:,:,:) / 216.0d0

      end do