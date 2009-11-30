c
c     ----------------------------------------------
c     CASE 7 - quarks at (11,0), (0,-6) and (0, 6)
c     ----------------------------------------------
c
c     base
c
      tl1r = cshift( bxlr, dim=ydir, shift= 6)
      tl1i = cshift( bxli, dim=ydir, shift= 6)

      sbxlr = bxlr
      sbxli = bxli

      bxlr = 0.0d0
      bxli = 0.0d0

      do ic = 1, nc
         do jc = 1, nc
            do kc = 1, nc

               bxlr(:,:,:,:,ic,jc) = bxlr(:,:,:,:,ic,jc) +
     &              sbxlr(:,:,:,:,ic,kc)*tl1r(:,:,:,:,kc,jc) - sbxli(:,:,:,:,ic,kc)*tl1i(:,:,:,:,kc,jc)

               bxli(:,:,:,:,ic,jc) = bxli(:,:,:,:,ic,jc) +
     &              sbxlr(:,:,:,:,ic,kc)*tl1i(:,:,:,:,kc,jc) + sbxli(:,:,:,:,ic,kc)*tl1r(:,:,:,:,kc,jc)

            end do
         end do
      end do
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

      tlr = cshift( cshift( tl1r, dim=ydir, shift= 4), dim=xdir, shift=-7)
      tli = cshift( cshift( tl1i, dim=ydir, shift= 4), dim=xdir, shift=-7)

      bllr = cshift( cshift( ld3r, dim=ydir, shift= 6),dim=xdir, shift=-4)
      blli = cshift( cshift( ld3i, dim=ydir, shift= 6),dim=xdir, shift=-4)

      tl2r = 0.0d0
      tl2i = 0.0d0

      do ic = 1, nc
         do jc = 1, nc
            do kc = 1, nc

               tl2r(:,:,:,:,ic,jc) = tl2r(:,:,:,:,ic,jc) +
     &              tl1r(:,:,:,:,ic,kc)*bllr(:,:,:,:,kc,jc) - tl1i(:,:,:,:,ic,kc)*blli(:,:,:,:,kc,jc)

               tl2i(:,:,:,:,ic,jc) = tl2i(:,:,:,:,ic,jc) +
     &              tl1r(:,:,:,:,ic,kc)*blli(:,:,:,:,kc,jc) + tl1i(:,:,:,:,ic,kc)*bllr(:,:,:,:,kc,jc)

            end do
         end do
      end do

      bllr = 0.d0
      blli = 0.d0

      do ic = 1, nc
         do jc = 1, nc
            do kc = 1, nc

               bllr(:,:,:,:,ic,jc) = bllr(:,:,:,:,ic,jc) +
     &              tl2r(:,:,:,:,ic,kc)*tlr(:,:,:,:,kc,jc) - tl2i(:,:,:,:,ic,kc)*tli(:,:,:,:,kc,jc)

               blli(:,:,:,:,ic,jc) = blli(:,:,:,:,ic,jc) +
     &              tl2r(:,:,:,:,ic,kc)*tli(:,:,:,:,kc,jc) + tl2i(:,:,:,:,ic,kc)*tlr(:,:,:,:,kc,jc)

            end do
         end do
      end do

      bllr = cshift( cshift( bllr, dim=ydir, shift= 6),dim=xdir, shift=11)
      blli = cshift( cshift( blli, dim=ydir, shift= 6),dim=xdir, shift=11)
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
     &              ld2r(:,:,:,:,ic,kc)*tlr(:,:,:,:,kc,jc) - rd2i(:,:,:,:,ic,kc)*tli(:,:,:,:,kc,jc)

               tl1i(:,:,:,:,ic,jc) = tl1i(:,:,:,:,ic,jc) +
     &              ld2r(:,:,:,:,ic,kc)*tli(:,:,:,:,kc,jc) + rd2i(:,:,:,:,ic,kc)*tlr(:,:,:,:,kc,jc)

            end do
         end do
      end do

      tlr = cshift( cshift( tl1r, dim=ydir, shift=-4), dim=xdir, shift=-7)
      tli = cshift( cshift( tl1i, dim=ydir, shift=-4), dim=xdir, shift=-7)

      brlr = cshift( cshift( rd3r, dim=ydir, shift=-6),dim=xdir, shift=-4)
      brli = cshift( cshift( rd3i, dim=ydir, shift=-6),dim=xdir, shift=-4)

      tl2r = 0.0d0
      tl2i = 0.0d0

      do ic = 1, nc
         do jc = 1, nc
            do kc = 1, nc

               tl2r(:,:,:,:,ic,jc) = tl2r(:,:,:,:,ic,jc) +
     &              tl1r(:,:,:,:,ic,kc)*brlr(:,:,:,:,kc,jc) - tl1i(:,:,:,:,ic,kc)*brli(:,:,:,:,kc,jc)

               tl2i(:,:,:,:,ic,jc) = tl2i(:,:,:,:,ic,jc) +
     &              tl1r(:,:,:,:,ic,kc)*brli(:,:,:,:,kc,jc) + tl1i(:,:,:,:,ic,kc)*brlr(:,:,:,:,kc,jc)

            end do
         end do
      end do

      brlr = 0.d0
      brli = 0.d0

      do ic = 1, nc
         do jc = 1, nc
            do kc = 1, nc

               brlr(:,:,:,:,ic,jc) = brlr(:,:,:,:,ic,jc) +
     &              tl2r(:,:,:,:,ic,kc)*tlr(:,:,:,:,kc,jc) - tl2i(:,:,:,:,ic,kc)*tli(:,:,:,:,kc,jc)

               brli(:,:,:,:,ic,jc) = brli(:,:,:,:,ic,jc) +
     &              tl2r(:,:,:,:,ic,kc)*tli(:,:,:,:,kc,jc) + tl2i(:,:,:,:,ic,kc)*tlr(:,:,:,:,kc,jc)

            end do
         end do
      end do

      brlr = cshift( cshift( brlr, dim=ydir, shift= 6),dim=xdir, shift=11)
      brli = cshift( cshift( brli, dim=ydir, shift= 6),dim=xdir, shift=11)
c
c     time links shifting
c
      do it = 1, nL(that)

         tlr(:,:,:,:,:,:) = tltr(:,:,:,:,:,:,it)
         tli(:,:,:,:,:,:) = tlti(:,:,:,:,:,:,it)

         tl1r = cshift( tlr, dim=ydir,shift= 12)
         tl1i = cshift( tli, dim=ydir,shift= 12)

         tl2r = cshift( cshift( tlr, dim=ydir,shift= 6),dim=xdir,shift= 11)
         tl2i = cshift( cshift( tli, dim=ydir,shift= 6),dim=xdir,shift= 11)

         include 'res.f'

         W(:,:,:,:,5,it) = Res(:,:,:,:) / 216.0d0

      end do