c
c     ----------------------------------------------
c     CASE 4 - quarks at ( 5,0), (0,-3) and (0, 3)
c     ----------------------------------------------
c     
c     top leg ( 5,0)
c
         call product(ur,ui,xdir,bxlr,bxli,5)

         ustr = cshift(uptr,dim=xdir,shift=5) 
         usti = cshift(upti,dim=xdir,shift=5)
         
         include 'TYfiles/topleg.f'
c
c     left leg, i.e., quark at (0, 3)
c
         call product(ur,ui,ydir,bllr,blli,3)

         ustr = cshift(uptr,dim=ydir,shift= 3) 
         usti = cshift(upti,dim=ydir,shift= 3)
         
         include 'TYfiles/leftleg.f'
c     
c     right leg, i.e., quark at (0,-3)
c
         brlr = cshift(bllr,dim=ydir,shift=-3) 
         brli = cshift(blli,dim=ydir,shift=-3) 

         ustr = cshift(uptr,dim=ydir,shift=-3) 
         usti = cshift(upti,dim=ydir,shift=-3)
         
         include 'TYfiles/rightleg.f'

         include 'TYfiles/res.f'
