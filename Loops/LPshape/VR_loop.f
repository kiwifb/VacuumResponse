      do Tee = 1, nL(that)

         call product(utr,uti,that,uptr,upti,Tee)

         call ly_loops(xhat,yhat,utr,uti,uptr,upti,Tee,WLxy(:,:,:,:,:,:,Tee))
         call ly_loops(xhat,zhat,utr,uti,uptr,upti,Tee,WLxz(:,:,:,:,:,:,Tee))

         do ilp = 1, nloop

            WLavg(ilp,Tee) = (Sum( WLxy(:,:,:,:,ilp,:,Tee)) +
     &                        Sum( WLxz(:,:,:,:,ilp,:,Tee))) / (8*volume)

         end do

      end do
