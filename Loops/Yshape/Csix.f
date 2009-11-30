c
c     ----------------------------------------------
c     CASE 6 - quarks at ( 5,0), (-4,-5) and (-4, 5)
c     ----------------------------------------------
c
c     top leg ( 5,0)
c
         call product(ur,ui,xdir,bxlr,bxli,5)

         ustr = cshift(uptr,dim=xdir,shift=5)
         usti = cshift(upti,dim=xdir,shift=5)

         include 'topleg.f'
c
c     left leg, i.e., quark at (-4, 5)
c
c     to (-1, 1)
c
         bllr = up1x1r
         blli = up1x1i
c
c     to (-3, 4)
c
         call multboxes(xdir,ydir,-1,1,bllr,blli,up2x3r,up2x3i)
c
c     to (-4, 5)
c
         call multboxes(xdir,ydir,-3,4,bllr,blli,up1x1r,up1x1i)

         ustr = cshift( cshift(uptr, dim=xdir,shift=-4), dim=ydir,shift= 5)
         usti = cshift( cshift(upti, dim=xdir,shift=-4), dim=ydir,shift= 5)

         include 'leftleg.f'
c
c     right leg, i.e., quark at (-4,-5)
c
c     to (-1,-1)
c
         brlr = do1x1r
         brli = do1x1i
c
c     to (-3,-4)
c
         call multboxes(xdir,ydir,-1,-1,brlr,brli,do2x3r,do2x3i)
c
c     to (-4,-5)
c
         call multboxes(xdir,ydir,-3,-4,brlr,brli,do1x1r,do1x1i)

         ustr = cshift( cshift(uptr, dim=xdir,shift=-4), dim=ydir,shift=-5)
         usti = cshift( cshift(upti, dim=xdir,shift=-4), dim=ydir,shift=-5)

         include 'rightleg.f'

         include 'res.f'
