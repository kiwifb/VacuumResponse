c
c     ----------------------------------------------
c     CASE 2 - quarks at ( 2,0), (-1,-2) and (-1, 2)
c     ----------------------------------------------
c
c     top leg ( 2,0)
c
         call product(ur,ui,xdir,bxlr,bxli,2)

         ustr = cshift(uptr,dim=xdir,shift=2) 
         usti = cshift(upti,dim=xdir,shift=2)
         
         include 'Yfiles/topleg.f'
c
c     left leg, i.e., quark at (-1, 2)
c
         bllr = ld2r
         blli = ld2i

         ustr = cshift(cshift(uptr,dim=xdir,shift=-1),dim=ydir,shift= 2) 
         usti = cshift(cshift(upti,dim=xdir,shift=-1),dim=ydir,shift= 2)
         
         include 'Yfiles/leftleg.f'
c
c     right leg, i.e., quark at (-1,-2)
c
         brlr = rd2r
         brli = rd2i

         ustr = cshift(cshift(uptr,dim=xdir,shift=-1),dim=ydir,shift=-2) 
         usti = cshift(cshift(upti,dim=xdir,shift=-1),dim=ydir,shift=-2)
         
         include 'Yfiles/rightleg.f'

         include 'Yfiles/res.f'
