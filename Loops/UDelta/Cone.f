c
c     ----------------------------------------------
c     CASE 1 - quarks at ( 2,0), (0,-1) and (0, 1)
c     ----------------------------------------------
c
c     base (bx)
c
      call product(ur,ui,ydir,cxlr,cxli,1)

      dxlr = cshift( cxlr, dim=ydir,shift=-1)
      dxli = cshift( cxli, dim=ydir,shift=-1)

      include 'basex.f'

      call product(ur,ui,xdir,cylr,cyli,1)
      call product(ur,ui,xdir,cylr,cyli,2)
c
c     top diagonal (bl)
c
      dylr = cshift( cylr, dim=ydir,shift= 1)
      dyli = cshift( cyli, dim=ydir,shift= 1)

      tmp1r = cshift( cxlr, dim=xdir,shift= 2)
      tmp1i = cshift( cxli, dim=xdir,shift= 2)

      include 'baseup.f'
c
c     bottom diagonal (br)
c
      dylr = cshift( cylr, dim=ydir,shift=-1)
      dyli = cshift( cyli, dim=ydir,shift=-1)

      tmp1r = cshift( dxlr, dim=xdir,shift= 2)
      tmp1i = cshift( dxli, dim=xdir,shift= 2)

      include 'basedo.f'
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
