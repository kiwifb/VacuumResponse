c
c     ----------------------------------------------
c     CASE 2 - quarks at ( 2,0), (-1,-2) and (-1, 2)
c     ----------------------------------------------
c
c     top leg ( 2,0)
c
         call product(ur,ui,xdir,bxlr,bxli,2)

         ustr = cshift(uptr,dim=xdir,shift=2)
         usti = cshift(upti,dim=xdir,shift=2)

         include 'topleg.f'
c
c     left/up leg, i.e., quark at (-1, 2)
c
         bllr = up1x2r
         blli = up1x2i

         ustr = cshift(cshift(uptr,dim=xdir,shift=-1),dim=ydir,shift= 2)
         usti = cshift(cshift(upti,dim=xdir,shift=-1),dim=ydir,shift= 2)

         include 'leftleg.f'
c
c     right/down leg, i.e., quark at (-1,-2)
c
         brlr = do1x2r
         brli = do1x2i

         ustr = cshift(cshift(uptr,dim=xdir,shift=-1),dim=ydir,shift=-2)
         usti = cshift(cshift(upti,dim=xdir,shift=-1),dim=ydir,shift=-2)

         include 'rightleg.f'

         include 'res.f'
