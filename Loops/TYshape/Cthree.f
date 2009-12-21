c
c     ----------------------------------------------
c     CASE 3 - quarks at ( 4,0), (0,-2) and (0, 2)
c     ----------------------------------------------
c
c     top leg ( 3,0)
c
         call product(ur,ui,xdir,bxlr,bxli,4)

         ustr = cshift(uptr,dim=xdir,shift=4)
         usti = cshift(upti,dim=xdir,shift=4)

         include 'topleg.f'
c
c     left leg, i.e., quark at (0, 2),
c     no change from previous case
c
c     right leg, i.e., quark at (0,-2)
c     no change from previous case
c
         include 'res.f'
