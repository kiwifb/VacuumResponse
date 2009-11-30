c
c     ----------------------------------------------
c     CASE 3 - quarks at ( 3,0), (-1,-2) and (-1, 2)
c     ----------------------------------------------
c     
c     top leg ( 3,0)
c     
         call product(ur,ui,xdir,bxlr,bxli,3)

         ustr = cshift(uptr,dim=xdir,shift=3) 
         usti = cshift(upti,dim=xdir,shift=3)
         
         include 'Yfiles/topleg.f'
c
c     left leg, i.e., quark at (-1, 2),
c     no change from previous case
c
c     right leg, i.e., quark at (-1,-2)
c     no change from previous case
c
         include 'Yfiles/res.f'
