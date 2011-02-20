      call Tetra_loops(utr,uti,WOP(:,:,:,:,:,:))

      do ilp = 1, nloop

        do it = 1, nL(that)

          WOP_avg(ilp,it)  = Sum(WOP(:,:,:,:,ilp,it) )  / volume

        end do ! it

      end do ! idq
