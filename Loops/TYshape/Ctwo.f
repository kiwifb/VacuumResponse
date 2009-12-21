c
c     ----------------------------------------------
c     CASE 2 - quarks at ( 3,0), (0,-2) and (0, 2)
c     ----------------------------------------------
c
c     top leg ( 3,0)
c
         call product(ur,ui,xdir,bxlr,bxli,3)

         ustr = cshift(uptr,dim=xdir,shift=3)
         usti = cshift(upti,dim=xdir,shift=3)

         include 'topleg.f'
c
c     left leg, i.e., quark at (0, 2)
c
         call product(ur,ui,ydir,bllr,blli,2)

         ustr = cshift(uptr,dim=ydir,shift= 2)
         usti = cshift(upti,dim=ydir,shift= 2)

         include 'leftleg.f'
c
c     right leg, i.e., quark at (0,-2)
c
         brlr = cshift(bllr,dim=ydir,shift=-2)
         brli = cshift(blli,dim=ydir,shift=-2)

         ustr = cshift(uptr,dim=ydir,shift=-2)
         usti = cshift(upti,dim=ydir,shift=-2)

         include 'rightleg.f'

         include 'res.f'
