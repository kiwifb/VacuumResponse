c
c     ----------------------------------------------
c     CASE 1 - quarks at ( 2,0), (0,-1) and (0, 1)
c     ----------------------------------------------
c
c     base (bx) 
c
      call product(ur,ui,ydir,bxlr,bxli,1)
      call product(ur,ui,ydir,bxlr,bxli,2)
c
c     top diagonal (bl)
c
      bllr = cshift( cshift( ld2r,dim=ydir,shift= 1),dim=xdir,shift=2)
      blli = cshift( cshift( ld2i,dim=ydir,shift= 1),dim=xdir,shift=2)
c
c     bottom diagonal (br)
c
      brlr = cshift( cshift( rd2r,dim=ydir,shift= 1),dim=xdir,shift=2)
      brli = cshift( cshift( rd2i,dim=ydir,shift= 1),dim=xdir,shift=2)
c
c     time links shifting
c
      do it = 1, nL(that)

         tlr(:,:,:,:,:,:) = tltr(:,:,:,:,:,:,it)
         tli(:,:,:,:,:,:) = tlti(:,:,:,:,:,:,it)

         tl1r = cshift( tlr, dim=ydir,shift= 2)
         tl1i = cshift( tli, dim=ydir,shift= 2)

         tl2r = cshift( cshift( tlr, dim=ydir,shift= 1),dim=xdir,shift= 2)
         tl2i = cshift( cshift( tli, dim=ydir,shift= 1),dim=xdir,shift= 2)

         include 'res.f'

         W(:,:,:,:,1,it) = Res(:,:,:,:) / 216.0d0

      end do