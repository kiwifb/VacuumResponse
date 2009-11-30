c
c     ----------------------------------------------
c     CASE 5 - quarks at ( 4,0), (-3,-4) and (-3, 4)
c     ----------------------------------------------
c
c     top leg ( 4,0)
c
         call product(ur,ui,xdir,bxlr,bxli,4)

         ustr = cshift(uptr,dim=xdir,shift=4)
         usti = cshift(upti,dim=xdir,shift=4)

         include 'topleg.f'
c
c     left/up leg, i.e., quark at (-3, 4)
c
c     to (-1, 1)
c
         bllr = up1x1r
         blli = up1x1i
c
c     to (-2, 3)
c
         call multboxes(xdir,ydir,-1,1,bllr,blli,up1x2r,up1x2i)
c
c     to (-3, 4)
c
         call multboxes(xdir,ydir,-2,3,bllr,blli,up1x1r,up1x1i)

         ustr = cshift(cshift(uptr,dim=xdir,shift=-3),dim=ydir,shift= 4)
         usti = cshift(cshift(upti,dim=xdir,shift=-3),dim=ydir,shift= 4)

         include 'leftleg.f'
c
c     right/down leg, i.e., quark at (-3,-4)
c
c     to (-1,-1)
c
         brlr = do1x1r
         brli = do1x1i
c
c     to (-2,-3)
c
         call multboxes(xdir,ydir,-1,-1,brlr,brli,do1x2r,do1x2i)
c
c     to (-3,-4)
c
         call multboxes(xdir,ydir,-2,-3,brlr,brli,do1x1r,do1x1i)

         ustr = cshift(cshift(uptr,dim=xdir,shift=-3),dim=ydir,shift=-4)
         usti = cshift(cshift(upti,dim=xdir,shift=-3),dim=ydir,shift=-4)

         include 'rightleg.f'

         include 'res.f'
