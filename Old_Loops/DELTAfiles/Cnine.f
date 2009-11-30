c
c     ----------------------------------------------
c     CASE 9 - quarks at (14,0), (0,-8) and (0, 8)
c     ----------------------------------------------
c     
c     top leg ( 14,0)
c     
         call product(ur,ui,xdir,bxlr,bxli,13)
         call product(ur,ui,xdir,bxlr,bxli,14)

         ustr = cshift(uptr,dim=xdir,shift=14)
         usti = cshift(upti,dim=xdir,shift=14)
         
         include 'topleg.f'
c
c     left leg, i.e., quark at (0, 8)
c
         call product(ur,ui,ydir,bllr,blli,8)

         ustr = cshift(uptr,dim=ydir,shift= 8) 
         usti = cshift(upti,dim=ydir,shift= 8)
         
         include 'leftleg.f'
c     
c     right leg, i.e., quark at (0,-8)
c
         brlr = cshift(bllr,dim=ydir,shift=-8) 
         brli = cshift(blli,dim=ydir,shift=-8) 
         
         ustr = cshift(uptr,dim=ydir,shift=-8) 
         usti = cshift(upti,dim=ydir,shift=-8)
         
         include 'rightleg.f'

         include 'res.f'         
