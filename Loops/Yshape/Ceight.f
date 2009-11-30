c
c     ----------------------------------------------
c     CASE 8 - quarks at ( 8,0), (-4,-7) and (-4, 7)
c     ----------------------------------------------
c
c     top leg ( 8,0)
c
         call product(ur,ui,xdir,bxlr,bxli,8)

         ustr = cshift(uptr,dim=xdir,shift=8)
         usti = cshift(upti,dim=xdir,shift=8)

         include 'topleg.f'
c
c     left leg, i.e., quark at (-4, 7)
c
c     to (-1, 2)
c
         bllr = up1x2r
         blli = up1x2i
c
c     to (-3, 5)
c
         call multboxes(xdir,ydir,-1,2,bllr,blli,up2x3r,up2x3i)
c
c     to (-4, 7)
c
         call multboxes(xdir,ydir,-3,5,bllr,blli,up1x2r,up1x2i)

         ustr = cshift( cshift(uptr, dim=xdir,shift=-4), dim=ydir,shift= 7)
         usti = cshift( cshift(upti, dim=xdir,shift=-4), dim=ydir,shift= 7)

         include 'leftleg.f'
c
c     right leg, i.e., quark at (-4,-7)
c
c     to (-1,-2)
c
         brlr = do1x2r
         brli = do1x2i
c
c     to (-3,-5)
c
         call multboxes(xdir,ydir,-1,-2,brlr,brli,do2x3r,do2x3i)
c
c     to (-4,-7)
c
         call multboxes(xdir,ydir,-3,-5,brlr,brli,do1x2r,do1x2i)

         ustr = cshift( cshift(uptr, dim=xdir,shift=-4), dim=ydir,shift=-7)
         usti = cshift( cshift(upti, dim=xdir,shift=-4), dim=ydir,shift=-7)

         include 'rightleg.f'

         include 'res.f'
