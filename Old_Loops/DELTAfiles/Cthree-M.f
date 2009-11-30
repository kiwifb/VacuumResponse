c
c     ----------------------------------------------
c     CASE 3 - quarks at ( 4,0), (0,-2) and (0, 2)
c     ----------------------------------------------
c
c     base unchanged
c
c     top diagonal
c
      tlr = cshift( cshift( ld2r, dim=ydir, shift= 1), dim=xdir, shift=-2)
      tli = cshift( cshift( ld2i, dim=ydir, shift= 1), dim=xdir, shift=-2)

      tl1r = 0.0d0
      tl1i = 0.0d0

      do ic = 1, nc
         do jc = 1, nc
            do kc = 1, nc

               tl1r(:,:,:,:,ic,jc) = tl1r(:,:,:,:,ic,jc) +
     &              ld2r(:,:,:,:,ic,kc)*tlr(:,:,:,:,kc,jc) - ld2i(:,:,:,:,ic,kc)*tli(:,:,:,:,kc,jc)

               tl1i(:,:,:,:,ic,jc) = tl1i(:,:,:,:,ic,jc) +
     &              ld2r(:,:,:,:,ic,kc)*tli(:,:,:,:,kc,jc) + ld2i(:,:,:,:,ic,kc)*tlr(:,:,:,:,kc,jc)

            end do
         end do
      end do

      bllr = cshift( cshift( tl1r, dim=ydir, shift= 2),dim=xdir, shift= 4)
      blli = cshift( cshift( tl1i, dim=ydir, shift= 2),dim=xdir, shift= 4)
c
c     bottom diagonal
c
      tlr = cshift( cshift( rd2r, dim=ydir, shift=-1), dim=xdir, shift=-2)
      tli = cshift( cshift( rd2i, dim=ydir, shift=-1), dim=xdir, shift=-2)

      tl1r = 0.0d0
      tl1i = 0.0d0

      do ic = 1, nc
         do jc = 1, nc
            do kc = 1, nc

               tl1r(:,:,:,:,ic,jc) = tl1r(:,:,:,:,ic,jc) +
     &              rd2r(:,:,:,:,ic,kc)*tlr(:,:,:,:,kc,jc) - rd2i(:,:,:,:,ic,kc)*tli(:,:,:,:,kc,jc)

               tl1i(:,:,:,:,ic,jc) = tl1i(:,:,:,:,ic,jc) +
     &              rd2r(:,:,:,:,ic,kc)*tli(:,:,:,:,kc,jc) + rd2i(:,:,:,:,ic,kc)*tlr(:,:,:,:,kc,jc)

            end do
         end do
      end do

      brlr = cshift( cshift( tl1r, dim=ydir, shift= 2),dim=xdir, shift= 4)
      brli = cshift( cshift( tl1i, dim=ydir, shift= 2),dim=xdir, shift= 4)
c
c     time links shifting
c
      do it = 1, nL(that)

         tlr(:,:,:,:,:,:) = tltr(:,:,:,:,:,:,it)
         tli(:,:,:,:,:,:) = tlti(:,:,:,:,:,:,it)

         tl1r = cshift( tlr, dim=ydir,shift= 4)
         tl1i = cshift( tli, dim=ydir,shift= 4)

         tl2r = cshift( cshift( tlr, dim=ydir,shift= 2),dim=xdir,shift= 4)
         tl2i = cshift( cshift( tli, dim=ydir,shift= 2),dim=xdir,shift= 4)

         include 'res.f'

         Wp(:,:,:,:,3,it) = Res(:,:,:,:) / 216.0d0

         call delta_mirror(xdir,ydir,bxlr,bxli,bllr,blli,brlr,brli,tlr,tli,tl1r,tl1i,Res,it,4,2)

         Wm(:,:,:,:,3,it) = Res(:,:,:,:) / 216.0d0

      end do