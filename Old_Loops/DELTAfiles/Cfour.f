c
c     ----------------------------------------------
c     CASE 4 - quarks at ( 5,0), (0,-3) and (0, 3)
c     ----------------------------------------------
c
c     base
c
      call product(ur,ui,ydir,bxlr,bxli,5)
      call product(ur,ui,ydir,bxlr,bxli,6)
c
c     top diagonal
c
      tlr = cshift( cshift( ld1r, dim=ydir, shift= 1), dim=xdir, shift=-2)
      tli = cshift( cshift( ld1i, dim=ydir, shift= 1), dim=xdir, shift=-2)

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

      tlr = cshift( cshift( ld2r, dim=ydir, shift= 2), dim=xdir, shift=-3)
      tli = cshift( cshift( ld2i, dim=ydir, shift= 2), dim=xdir, shift=-3)

      tl2r = 0.0d0
      tl2i = 0.0d0

      do ic = 1, nc
         do jc = 1, nc
            do kc = 1, nc

               tl2r(:,:,:,:,ic,jc) = tl2r(:,:,:,:,ic,jc) +
     &              tl1r(:,:,:,:,ic,kc)*tlr(:,:,:,:,kc,jc) - tl1i(:,:,:,:,ic,kc)*tli(:,:,:,:,kc,jc)

               tl2i(:,:,:,:,ic,jc) = tl2i(:,:,:,:,ic,jc) +
     &              tl1r(:,:,:,:,ic,kc)*tli(:,:,:,:,kc,jc) + tl1i(:,:,:,:,ic,kc)*tlr(:,:,:,:,kc,jc)

            end do
         end do
      end do

      bllr = cshift( cshift( tl2r, dim=ydir, shift= 3),dim=xdir, shift= 5)
      blli = cshift( cshift( tl2i, dim=ydir, shift= 3),dim=xdir, shift= 5)
c
c     bottom diagonal
c
      tlr = cshift( cshift( rd1r, dim=ydir, shift=-1), dim=xdir, shift=-2)
      tli = cshift( cshift( rd1i, dim=ydir, shift=-1), dim=xdir, shift=-2)

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

      tlr = cshift( cshift( rd2r, dim=ydir, shift=-2), dim=xdir, shift=-3)
      tli = cshift( cshift( rd2i, dim=ydir, shift=-2), dim=xdir, shift=-3)

      tl2r = 0.0d0
      tl2i = 0.0d0

      do ic = 1, nc
         do jc = 1, nc
            do kc = 1, nc

               tl2r(:,:,:,:,ic,jc) = tl2r(:,:,:,:,ic,jc) +
     &              tl1r(:,:,:,:,ic,kc)*tlr(:,:,:,:,kc,jc) - tl1i(:,:,:,:,ic,kc)*tli(:,:,:,:,kc,jc)

               tl2i(:,:,:,:,ic,jc) = tl2i(:,:,:,:,ic,jc) +
     &              tl1r(:,:,:,:,ic,kc)*tli(:,:,:,:,kc,jc) + tl1i(:,:,:,:,ic,kc)*tlr(:,:,:,:,kc,jc)

            end do
         end do
      end do

      brlr = cshift( cshift( tl2r, dim=ydir, shift= 3),dim=xdir, shift= 5)
      brli = cshift( cshift( tl2i, dim=ydir, shift= 3),dim=xdir, shift= 5)
c
c     time links shifting
c
      do it = 1, nL(that)

         tlr(:,:,:,:,:,:) = tltr(:,:,:,:,:,:,it)
         tli(:,:,:,:,:,:) = tlti(:,:,:,:,:,:,it)

         tl1r = cshift( tlr, dim=ydir,shift= 6)
         tl1i = cshift( tli, dim=ydir,shift= 6)

         tl2r = cshift( cshift( tlr, dim=ydir,shift= 3),dim=xdir,shift= 5)
         tl2i = cshift( cshift( tli, dim=ydir,shift= 3),dim=xdir,shift= 5)

         include 'res.f'

         W(:,:,:,:,4,it) = Res(:,:,:,:) / 216.0d0

      end do