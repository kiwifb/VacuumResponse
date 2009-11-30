c
c     ----------------------------------------------
c     CASE 8 - quarks at (12,0), (0,-7) and (0, 7)
c     ----------------------------------------------
c
c     top leg (12,0)
c
         call product(ur,ui,xdir,bxlr,bxli,12)

         ustr = cshift(uptr,dim=xdir,shift=12)
         usti = cshift(upti,dim=xdir,shift=12)

         include 'topleg.f'
c
c     left leg, i.e., quark at (0, 7)
c
         call product(ur,ui,ydir,bllr,blli,7)

         ustr = cshift(uptr,dim=ydir,shift= 7) 
         usti = cshift(upti,dim=ydir,shift= 7)

         include 'leftleg.f'
c
c     right leg, i.e., quark at (0,-7)
c
         brlr = cshift(bllr,dim=ydir,shift=-7) 
         brli = cshift(blli,dim=ydir,shift=-7) 

         ustr = cshift(uptr,dim=ydir,shift=-7) 
         usti = cshift(upti,dim=ydir,shift=-7)

         include 'rightleg.f'

         include 'res.f'
