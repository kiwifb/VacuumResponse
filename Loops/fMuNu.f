c
c     fMuNu calculates the field strength tensor via lotsaloops
c     Made into a module by F. Bissey Dec. 2003
c
      MODULE L_fMuNu

      CONTAINS

      subroutine fMuNu(ur,ui,Fr,Fi,uzero,Qpaths)

      USE GS_LATTICESIZE
      USE L_lotsaLoops

      implicit none

c     global variables

      integer,parameter                                       :: nc=3                   !no. colours
      integer,parameter                                       :: mu=4                   !no. dimensions

      integer,parameter                                       :: nf=6

      double precision,dimension(nx,ny,nz,nt,mu,nc,nc)        :: ur,ui
!HPF$ DISTRIBUTE ur(*,*,BLOCK,BLOCK,*,*,*)
!HPF$ DISTRIBUTE ui(*,*,BLOCK,BLOCK,*,*,*)

      double precision,dimension(nx,ny,nz,nt,nc,nc,nf)        :: Fr,Fi                  !Fmunu
!HPF$ DISTRIBUTE Fr(*,*,BLOCK,BLOCK,*,*,*)
!HPF$ DISTRIBUTE Fi(*,*,BLOCK,BLOCK,*,*,*)

      double precision                                        :: uzero
      integer                                                 :: Qpaths

c     local variables

      double precision,dimension(nx,ny,nz,nt,nc,nc)           :: Cr,Ci                  !Clover Term
!HPF$ DISTRIBUTE Cr(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE Ci(*,*,BLOCK,BLOCK,*,*)
     
      double precision,dimension(nx,ny,nz,nt)                 :: TrFr                   !Trace Variable
!HPF$ DISTRIBUTE TrFr(*,*,BLOCK,BLOCK)

      double precision,dimension(nx,ny,nz,nt,nc,nc)           :: stapler,staplei
!HPF$ DISTRIBUTE stapler(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE staplei(*,*,BLOCK,BLOCK,*,*)

      double precision,dimension(nx,ny,nz,nt,mu,nc,nc)           :: fsr,fsi
!HPF$ DISTRIBUTE fsr(*,*,BLOCK,BLOCK,*,*,*)
!HPF$ DISTRIBUTE fsi(*,*,BLOCK,BLOCK,*,*,*)

      integer                                                 :: ic,jc,kc,muhat,nuhat   !loop variables
      integer                                                 :: index

c     start of the execution commands.
c
c     Calculation of F_mu_nu
c
      Fr = 0.0d0
      Fi = 0.0d0

      do muhat=2,mu
         do nuhat=1,muhat-1
c
c     ReZero the clover term
c
            Cr = 0.0d0
            Ci = 0.0d0

c     Calculate the index that represents the unique combination of
c     muhat and nuhat.
c     This is done to save storage.
c     F(x,y,z,t,a,b,mu,nu) contains 16 elements but we can exploit the
c     traceless antisymmetric properties of F and only store 6
c     elements for each F.
c     ie
c     F(x,y,z,t,a,b,2,1)->F(x,a,b,2)
c     F(x,y,z,t,a,b,3,1)->F(x,a,b,3)
c     F(x,y,z,t,a,b,3,2)->F(x,a,b,4)
c     F(x,y,z,t,a,b,4,1)->F(x,a,b,1)
c     F(x,y,z,t,a,b,4,2)->F(x,a,b,5)
c     F(x,y,z,t,a,b,4,3)->F(x,a,b,6)
c     The unique muhat/nuhat index combination is rolled into the
c     single index by the following relation

            index = ( muhat + nuhat ) - 3 * ( ( muhat-nuhat ) / 3 ) - 1

c     Calculation of link products for upper-right plaquette, #1

c     t1 = U(x,mu)*U(x+mu,nu)
c
            stapler = 0.0d0
            staplei = 0.0d0

            call LotsaLoops(ur,ui,stapler,staplei,muhat,1,Qpaths,uzero,nuhat,.true.)

            do ic=1,nc
               do jc=1,nc
                  do kc=1,nc

                     Cr(:,:,:,:,ic,jc) = Cr(:,:,:,:,ic,jc)     +
     &                    ( ur(:,:,:,:,muhat,ic,kc) * stapler(:,:,:,:,kc,jc) 
     &                   -  ui(:,:,:,:,muhat,ic,kc) * staplei(:,:,:,:,kc,jc) ) / 4.0d0
                     Ci(:,:,:,:,ic,jc) = Ci(:,:,:,:,ic,jc)     +
     &                    ( ur(:,:,:,:,muhat,ic,kc) * staplei(:,:,:,:,kc,jc)
     &                   +  ui(:,:,:,:,muhat,ic,kc) * stapler(:,:,:,:,kc,jc) ) / 4.0d0

                  end do
               end do
            end do
c
c     Calculation of link products for upper-left plaquette, #2
c
            fsr = cshift(ur(:,:,:,:,:,:,:), dim = muhat, shift = -1)
            fsi = cshift(ui(:,:,:,:,:,:,:), dim = muhat, shift = -1)

            stapler = 0.0d0
            staplei = 0.0d0

            call LotsaLoops(fsr,fsi,stapler,staplei,muhat,2,Qpaths,uzero,nuhat,.true.)

            do ic=1,nc
               do jc=1,nc
                  do kc=1,nc

                     Cr(:,:,:,:,ic,jc) = Cr(:,:,:,:,ic,jc)     +
     &                    ( stapler(:,:,:,:,ic,kc) * fsr(:,:,:,:,muhat,kc,jc)
     &                   -  staplei(:,:,:,:,ic,kc) * fsi(:,:,:,:,muhat,kc,jc) ) / 4.0d0
                     Ci(:,:,:,:,ic,jc) = Ci(:,:,:,:,ic,jc)     +
     &                    ( stapler(:,:,:,:,ic,kc) * fsi(:,:,:,:,muhat,kc,jc)
     &                   +  staplei(:,:,:,:,ic,kc) * fsr(:,:,:,:,muhat,kc,jc) ) / 4.0d0

                  end do
               end do
            end do
c
c     Calculation of link products for lower-left plaquette, #3
c
            stapler = 0.0d0
            staplei = 0.0d0

            call LotsaLoops(fsr,fsi,stapler,staplei,muhat,3,Qpaths,uzero,nuhat,.true.)

            do ic=1,nc
               do jc=1,nc
                  do kc=1,nc

                     Cr(:,:,:,:,ic,jc) = Cr(:,:,:,:,ic,jc)     +
     &                    ( fsr(:,:,:,:,muhat,kc,ic) * stapler(:,:,:,:,jc,kc)
     &                    - fsi(:,:,:,:,muhat,kc,ic) * staplei(:,:,:,:,jc,kc) ) / 4.0d0
                     Ci(:,:,:,:,ic,jc) = Ci(:,:,:,:,ic,jc)     +
     &                  ( - fsr(:,:,:,:,muhat,kc,ic) * staplei(:,:,:,:,jc,kc)
     &                    - fsi(:,:,:,:,muhat,kc,ic) * stapler(:,:,:,:,jc,kc) ) / 4.0d0

                  end do
               end do
            end do
c
c     Calculation of link products for lower-right plaquette, #4
c
            stapler = 0.0d0
            staplei = 0.0d0

            call LotsaLoops(ur,ui,stapler,staplei,muhat,4,Qpaths,uzero,nuhat,.true.)

            do ic=1,nc
               do jc=1,nc
                  do kc=1,nc

                     Cr(:,:,:,:,ic,jc) = Cr(:,:,:,:,ic,jc)     +
     &                    ( stapler(:,:,:,:,kc,ic) * ur(:,:,:,:,muhat,jc,kc)
     &                    - staplei(:,:,:,:,kc,ic) * ui(:,:,:,:,muhat,jc,kc) ) / 4.0d0
                     Ci(:,:,:,:,ic,jc) = Ci(:,:,:,:,ic,jc)     +
     &                  ( - stapler(:,:,:,:,kc,ic) * ui(:,:,:,:,muhat,jc,kc)
     &                    - staplei(:,:,:,:,kc,ic) * ur(:,:,:,:,muhat,jc,kc) ) / 4.0d0

                  end do
               end do
            end do
c
c     Now lets implement the Hermitian and traceless aspects of
c     the Gell-Mann Matrices by subtracting the 
c     Hermitian conjugate of C which gives us a clover term of  C = 2i*g*F_mu_nu.  
c     Hence multiply C by -i/2 to get g*F_mu_nu. 
c
            do ic=1,nc
               do jc=1,nc

                  Fi(:,:,:,:,ic,jc,index) = -(Cr(:,:,:,:,ic,jc) - Cr(:,:,:,:,jc,ic))/2.0d0
                  Fr(:,:,:,:,ic,jc,index) =  (Ci(:,:,:,:,ic,jc) + Ci(:,:,:,:,jc,ic))/2.0d0

               end do
            end do
c
c     Fi is traceless, Fr is not, but the Gell-Mann Matrices are. 
c     Therefore subtract 1/3 the trace from each diagonal element of Fr.
c
            TrFr = 0.0d0

            do ic=1,nc

               TrFr(:,:,:,:) = TrFr(:,:,:,:) + Fr(:,:,:,:,ic,ic,index)

            end do

            do ic=1,nc

               Fr(:,:,:,:,ic,ic,index) = Fr(:,:,:,:,ic,ic,index) - TrFr(:,:,:,:)/3.0d0

            end do

         end do
      end do

      return

      end subroutine fMuNu

      END module L_fMuNu
