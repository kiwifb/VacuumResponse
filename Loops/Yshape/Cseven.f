c
c     ----------------------------------------------
c     CASE 7 - quarks at ( 7,0), (-4,-6) and (-4, 6)
c     ----------------------------------------------
c
c     top leg ( 7,0)
c
         call product(ur,ui,xdir,bxlr,bxli,6)
         call product(ur,ui,xdir,bxlr,bxli,7)

         ustr = cshift(uptr,dim=xdir,shift=7)
         usti = cshift(upti,dim=xdir,shift=7)

         include 'topleg.f'
c
c     left leg, i.e., quark at (-4, 6)
c
c     to (-2, 3)
c
         bllr = up2x3r
         blli = up2x3i
c
c     to (-4, 6)
c
         call multboxes(xdir,ydir,-2,3,bllr,blli,up2x3r,up2x3i)

         ustr = cshift( cshift(uptr, dim=xdir,shift=-4), dim=ydir,shift= 6)
         usti = cshift( cshift(upti, dim=xdir,shift=-4), dim=ydir,shift= 6)

         include 'leftleg.f'
c
c     right leg, i.e., quark at (-4,-6)
c
c     to (-2,-3)
c
         brlr = do2x3r
         brli = do2x3i
c
c     to (-4,-6)
c
         call multboxes(xdir,ydir,-2,-3,brlr,brli,do2x3r,do2x3i)

         ustr = cshift( cshift(uptr, dim=xdir,shift=-4), dim=ydir,shift=-6)
         usti = cshift( cshift(upti, dim=xdir,shift=-4), dim=ydir,shift=-6)

         include 'rightleg.f'

         include 'res.f'
