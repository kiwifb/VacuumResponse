c
c     ----------------------------------------------
c     CASE 1 - quarks at ( 1,0), (-1,-1) and (-1, 1)
c     ----------------------------------------------
c
c     top leg ( 1,0)
c
         call product(ur,ui,xdir,bxlr,bxli,1)

         ustr = cshift(uptr,dim=xdir,shift=1) 
         usti = cshift(upti,dim=xdir,shift=1)

         include 'Yfiles/topleg.f'
c
c     left leg, i.e., quark at (-1, 1)
c
         bllr = ld1r
         blli = ld1i

         ustr = cshift(cshift(uptr,dim=xdir,shift=-1),dim=ydir,shift= 1) 
         usti = cshift(cshift(upti,dim=xdir,shift=-1),dim=ydir,shift= 1)

         include 'Yfiles/leftleg.f'
c
c     right leg, i.e., quark at (-1,-1)
c
         brlr = rd1r
         brli = rd1i

         ustr = cshift(cshift(uptr,dim=xdir,shift=-1),dim=ydir,shift=-1) 
         usti = cshift(cshift(upti,dim=xdir,shift=-1),dim=ydir,shift=-1)

         include 'Yfiles/rightleg.f'

         include 'Yfiles/res.f'
