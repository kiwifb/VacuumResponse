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
c     bottom diagonal (br)
c
      brlr = cshift( do2x1r, dim=xdir,shift=2)
      brli = cshift( do2x1i, dim=xdir,shift=2)
c
c     loops creation
c
      do it= 1, nL(that)
         call system_clock(timea)

c
c        shifting the time link into place
c
         t_1r(:,:,:,:,:,:) = cshift( tlr(:,:,:,:,:,:,it), dim=ydir,shift= 1)
         t_1i(:,:,:,:,:,:) = cshift( tli(:,:,:,:,:,:,it), dim=ydir,shift= 1)

         t_2r(:,:,:,:,:,:) = cshift( tlr(:,:,:,:,:,:,it), dim=xdir,shift= 2)
         t_2i(:,:,:,:,:,:) = cshift( tli(:,:,:,:,:,:,it), dim=xdir,shift= 2)

         t_3r(:,:,:,:,:,:) = cshift( tlr(:,:,:,:,:,:,it), dim=ydir,shift=-1)
         t_3i(:,:,:,:,:,:) = cshift( tli(:,:,:,:,:,:,it), dim=ydir,shift=-1)

         include 'res.f'

         W(:,:,:,:,1,it) = Res(:,:,:,:) / 6.0d0
         call system_clock(timeb)
         write(*,'(a,i2,f15.8,a)') 'Tee:',it,(real(timeb-timea)/real(count_rate)),' s'

      end do
