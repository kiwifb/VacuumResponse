c
c     ----------------------------------------------
c     CASE 4 - quarks at ( 3,0), (-2,-3) and (-2, 3)
c     ----------------------------------------------
c
c     top leg ( 3,0)
c     no change from previous case
c
c     left/up leg, i.e., quark at (-2, 3)
c
         bllr = up2x3r
         blli = up2x3i

         ustr = cshift(cshift(uptr,dim=xdir,shift=-2),dim=ydir,shift= 3)
         usti = cshift(cshift(upti,dim=xdir,shift=-2),dim=ydir,shift= 3)

         include 'leftleg.f'
c
c     right/down leg, i.e., quark at (-2,-3)
c
         brlr = do2x3r
         brli = do2x3i

         ustr = cshift(cshift(uptr,dim=xdir,shift=-2),dim=ydir,shift=-3)
         usti = cshift(cshift(upti,dim=xdir,shift=-2),dim=ydir,shift=-3)

         include 'rightleg.f'

         include 'res.f'
