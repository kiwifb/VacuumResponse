c
c     ----------------------------------------------
c     CASE 4 - quarks at ( 3,0), (-2,-3) and (-2, 3)
c     ----------------------------------------------
c     
c     top leg ( 3,0)
c     no change from previous case
c
c     left leg, i.e., quark at (-2, 3)
c
         bllr = ld3r
         blli = ld3i

         ustr = cshift(cshift(uptr,dim=xdir,shift=-2),dim=ydir,shift= 3) 
         usti = cshift(cshift(upti,dim=xdir,shift=-2),dim=ydir,shift= 3)
         
         include 'Yfiles/leftleg.f'
c     
c     right leg, i.e., quark at (-2,-3)
c
         brlr = rd3r
         brli = rd3i

         ustr = cshift(cshift(uptr,dim=xdir,shift=-2),dim=ydir,shift=-3) 
         usti = cshift(cshift(upti,dim=xdir,shift=-2),dim=ydir,shift=-3)
         
         include 'Yfiles/rightleg.f'

         include 'Yfiles/res.f'
