      call Planar_loops(xhat,yhat,utr,uti,WPLxy(:,:,:,:,:,:,:))
c      call Planar_loops(xhat,zhat,utr,uti,WPLxz(:,:,:,:,:,:,:))
c      call OPlan_loops(xhat,yhat,zhat,utr,uti,WOP(:,:,:,:,:,:,:))

      do idq = 1, dqsep
        do iqq = 1, ysep

          do it = 1, nL(that)

             WPLxy_avg(iqq,idq,it)  = Sum( WPLxy(:,:,:,:,iqq,idq,it) )  / volume
c             WPLxz_avg(iqq,idq,it)  = Sum( WPLxy(:,:,:,:,iqq,idq,it) )  / volume
c               WOP_avg(iqq,idq,it)  = Sum(    WOP(:,:,:,:,iqq,idq,it) )  / volume

          end do ! it

        end do ! iqq
      end do ! idq
