c
c     ----------------------------------------------
c     CASE 5 - quarks at ( 7,0), (0,-4) and (0, 4)
c     ----------------------------------------------
c     
c     top leg ( 7,0)
c     
         call product(ur,ui,xdir,bxlr,bxli,6)
         call product(ur,ui,xdir,bxlr,bxli,7)

         ustr = cshift(uptr,dim=xdir,shift=7)
         usti = cshift(upti,dim=xdir,shift=7)
         
         include 'topleg.f'
c
c     left leg, i.e., quark at (0, 4)
c
         call product(ur,ui,ydir,bllr,blli,4)
         
         ustr = cshift(uptr,dim=ydir,shift= 4) 
         usti = cshift(upti,dim=ydir,shift= 4)
         
         include 'leftleg.f'
c     
c     right leg, i.e., quark at (0,-4)
c
         brlr = cshift(bllr,dim=ydir,shift=-4) 
         brli = cshift(blli,dim=ydir,shift=-4) 

         ustr = cshift(uptr,dim=ydir,shift=-4) 
         usti = cshift(upti,dim=ydir,shift=-4)
         
         include 'rightleg.f'

         include 'res.f'         
