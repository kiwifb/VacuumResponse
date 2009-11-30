c
c     ----------------------------------------------
c     CASE 9 - quarks at ( 9,0), (-5,-8) and (-5, 8)
c     ----------------------------------------------
c
c     top leg ( 9,0)
c
         call product(ur,ui,xdir,bxlr,bxli,9)

         ustr = cshift(uptr,dim=xdir,shift=9)
         usti = cshift(upti,dim=xdir,shift=9)

         include 'topleg.f'
c
c     left leg, i.e., quark at (-5, 8)
c
c     to (-2, 3)
c
         bllr = up2x3r
         blli = up2x3i
c
c     to (-3, 5)
c
         call multboxes(xdir,ydir,-2,3,bllr,blli,up1x2r,up1x2i)
c
c     to (-5, 8)
c
         call multboxes(xdir,ydir,-3,5,bllr,blli,up2x3r,up2x3i)

         ustr = cshift( cshift(uptr, dim=xdir,shift=-5), dim=ydir,shift= 8)
         usti = cshift( cshift(upti, dim=xdir,shift=-5), dim=ydir,shift= 8)

         include 'leftleg.f'
c
c     right leg, i.e., quark at (-5,-8)
c
c     to (-2,-3)
c
         brlr = do2x3r
         brli = do2x3i
c
c     to (-3,-5)
c
         call multboxes(xdir,ydir,-2,-3,brlr,brli,do1x2r,do1x2i)
c
c     to (-5,-8)
c
         call multboxes(xdir,ydir,-3,-5,brlr,brli,do2x3r,do2x3i)

         ustr = cshift( cshift(uptr, dim=xdir,shift=-5), dim=ydir,shift=-8)
         usti = cshift( cshift(upti, dim=xdir,shift=-5), dim=ydir,shift=-8)

         include 'rightleg.f'

         include 'res.f'
