c
c     ----------------------------------------------
c     CASE 1 - quarks at ( 2,0), (0,-1) and (0, 1)
c     ----------------------------------------------
c
c     top leg ( 2,0)
c
         call product(ur,ui,xdir,bxlr,bxli,1)
         call product(ur,ui,xdir,bxlr,bxli,2)

         ustr = cshift(uptr,dim=xdir,shift=2)
         usti = cshift(upti,dim=xdir,shift=2)

         include 'topleg.f'
c
c     left leg, i.e., quark at (0, 1)
c
         call product(ur,ui,ydir,bllr,blli,1)

         ustr = cshift(uptr,dim=ydir,shift= 1)
         usti = cshift(upti,dim=ydir,shift= 1)

         include 'leftleg.f'
c
c     right leg, i.e., quark at (0,-1)
c
         brlr = cshift(bllr,dim=ydir,shift=-1)
         brli = cshift(blli,dim=ydir,shift=-1)

         ustr = cshift(uptr,dim=ydir,shift=-1)
         usti = cshift(upti,dim=ydir,shift=-1)

         include 'rightleg.f'

         include 'res.f'
