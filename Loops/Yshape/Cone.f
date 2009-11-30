c
c     ----------------------------------------------
c     CASE 1 - quarks at ( 1,0), (-1,-1) and (-1, 1)
c     ----------------------------------------------
c
c     top leg ( 1,0)
c
         call product(ur,ui,xdir,bxlr,bxli,1)

         ustr = cshift(uptr,dim=xdir,shift=1)
         usti = cshift(upti,dim=xdir,shift=1)

         include 'topleg.f'
c
c     left/up leg, i.e., quark at (-1, 1)
c
         bllr = up1x1r
         blli = up1x1i

         ustr = cshift(cshift(uptr,dim=xdir,shift=-1),dim=ydir,shift= 1)
         usti = cshift(cshift(upti,dim=xdir,shift=-1),dim=ydir,shift= 1)

         include 'leftleg.f'
c
c     right/down leg, i.e., quark at (-1,-1)
c
         brlr = do1x1r
         brli = do1x1i

         ustr = cshift(cshift(uptr,dim=xdir,shift=-1),dim=ydir,shift=-1)
         usti = cshift(cshift(upti,dim=xdir,shift=-1),dim=ydir,shift=-1)

         include 'rightleg.f'

         include 'res.f'
