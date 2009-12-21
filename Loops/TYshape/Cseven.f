c
c     ----------------------------------------------
c     CASE 7 - quarks at (11,0), (0,-6) and (0, 6)
c     ----------------------------------------------
c
c     top leg (11,0)
c
         call product(ur,ui,xdir,bxlr,bxli,10)
         call product(ur,ui,xdir,bxlr,bxli,11)

         ustr = cshift(uptr,dim=xdir,shift=11)
         usti = cshift(upti,dim=xdir,shift=11)

         include 'topleg.f'
c
c     left leg, i.e., quark at (0, 6)
c
         call product(ur,ui,ydir,bllr,blli,6)

         ustr = cshift(uptr,dim=ydir,shift= 6)
         usti = cshift(upti,dim=ydir,shift= 6)

         include 'leftleg.f'
c
c     right leg, i.e., quark at (0,-6)
c
         brlr = cshift(bllr,dim=ydir,shift=-6)
         brli = cshift(blli,dim=ydir,shift=-6)

         ustr = cshift(uptr,dim=ydir,shift=-6)
         usti = cshift(upti,dim=ydir,shift=-6)

         include 'rightleg.f'

         include 'res.f'
