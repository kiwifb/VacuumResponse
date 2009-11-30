      do ilp = 1, nloop

         actionLC( iy, iz, it,ilp,Tee,offset) = actionLC( iy, iz, it,ilp,Tee,offset) +
     &             sum( WLxy(:,:,:,:,ilp,1,Tee) * actionTz(:,:,:,:) )

!         topChgLC( iy, iz, it,ilp,Tee,offset) = topChgLC( iy, iz, it,ilp,Tee,offset) +
!     &             sum( WLxy(:,:,:,:,ilp,1,Tee) * topChgTz(:,:,:,:) )

         ict = modulo(nt - it , nt)

         actionLC( iy, iz,ict,ilp,Tee,offset) = actionLC( iy, iz,ict,ilp,Tee,offset) +
     &             sum( WLxy(:,:,:,:,ilp,4,Tee) * actionTz(:,:,:,:) )

!         topChgLC( iy, iz,ict,ilp,Tee,offset) = topChgLC( iy, iz,ict,ilp,Tee,offset) +
!     &            sum( WLxy(:,:,:,:,ilp,4,Tee) * topChgTz(:,:,:,:) )

         icz = modulo(nz - iz , nz)

         actionLC( iy,icz, it,ilp,Tee,offset) = actionLC( iy,icz, it,ilp,Tee,offset) +
     &            sum( WLxy(:,:,:,:,ilp,2,Tee) * actionTz(:,:,:,:) )

!         topChgLC( iy,icz, it,ilp,Tee,offset) = topChgLC( iy,icz, it,ilp,Tee,offset) +
!     &            sum( WLxy(:,:,:,:,ilp,2,Tee) * topChgTz(:,:,:,:) )

         actionLC( iy,icz,ict,ilp,Tee,offset) = actionLC( iy,icz,ict,ilp,Tee,offset) +
     &            sum( WLxy(:,:,:,:,ilp,3,Tee) * actionTz(:,:,:,:) )

!         topChgLC( iy,icz,ict,ilp,Tee,offset) = topChgLC( iy,icz,ict,ilp,Tee,offset) +
!     &            sum( WLxy(:,:,:,:,ilp,3,Tee) * topChgTz(:,:,:,:) )

         icz = modulo(     iy , nz)
         icy = modulo(ny - iz , ny)

         actionLC(icy,icz, it,ilp,Tee,offset) = actionLC(icy,icz, it,ilp,Tee,offset) +
     &            sum( WLxz(:,:,:,:,ilp,1,Tee) * actionTz(:,:,:,:) )

!         topChgLC(icy,icz, it,ilp,Tee,offset) = topChgLC(icy,icz, it,ilp,Tee,offset) +
!     &            sum( WLxz(:,:,:,:,ilp,1,Tee) * topChgTz(:,:,:,:) )

         ict  = modulo(nt - it , nt)

         actionLC(icy,icz,ict,ilp,Tee,offset) = actionLC(icy,icz,ict,ilp,Tee,offset) +
     &            sum( WLxz(:,:,:,:,ilp,4,Tee) * actionTz(:,:,:,:) )

!         topChgLC(icy,icz,ict,ilp,Tee,offset) = topChgLC(icy,icz,ict,ilp,Tee,offset) +
!     &            sum( WLxz(:,:,:,:,ilp,4,Tee) * topChgTz(:,:,:,:) )

         icz = modulo(nz - iy , nz)
         icy = modulo(     iz , ny)

         actionLC(icy,icz, it,ilp,Tee,offset) = actionLC(icy,icz, it,ilp,Tee,offset) +
     &            sum( WLxz(:,:,:,:,ilp,2,Tee) * actionTz(:,:,:,:) )

!         topChgLC(icy,icz, it,ilp,Tee,offset) = topChgLC(icy,icz, it,ilp,Tee,offset) +
!     &            sum( WLxz(:,:,:,:,ilp,2,Tee) * topChgTz(:,:,:,:) )

         actionLC(icy,icz,ict,ilp,Tee,offset) = actionLC(icy,icz,ict,ilp,Tee,offset) +
     &            sum( WLxz(:,:,:,:,ilp,3,Tee) * actionTz(:,:,:,:) )

!         topChgLC(icy,icz, it,ilp,Tee,offset) = topChgLC(icy,icz, it,ilp,Tee,offset) +
!     &            sum( WLxz(:,:,:,:,ilp,3,Tee) * topChgTz(:,:,:,:) )

      end do     ! ilp loop end
