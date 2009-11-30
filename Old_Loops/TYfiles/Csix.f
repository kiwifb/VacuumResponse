c
c     ----------------------------------------------
c     CASE 6 - quarks at ( 9,0), (0,-5) and (0, 5)
c     ----------------------------------------------
c     
c     top leg ( 9,0)
c     
         call product(ur,ui,xdir,bxlr,bxli,8)
         call product(ur,ui,xdir,bxlr,bxli,9)

         ustr = cshift(uptr,dim=xdir,shift=9)
         usti = cshift(upti,dim=xdir,shift=9)
         
         include 'TYfiles/topleg.f'
c
c     left leg, i.e., quark at (0, 5)
c
         call product(ur,ui,ydir,bllr,blli,5)

         ustr = cshift(uptr,dim=ydir,shift= 5) 
         usti = cshift(upti,dim=ydir,shift= 5)
         
         include 'TYfiles/leftleg.f'
c     
c     right leg, i.e., quark at (0,-5)
c
         brlr = cshift(bllr,dim=ydir,shift=-5) 
         brli = cshift(blli,dim=ydir,shift=-5) 
         
         ustr = cshift(uptr,dim=ydir,shift=-5) 
         usti = cshift(upti,dim=ydir,shift=-5)
         
         include 'TYfiles/rightleg.f'

         include 'TYfiles/res.f'
