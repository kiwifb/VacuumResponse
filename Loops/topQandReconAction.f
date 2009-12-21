c     cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     TopQandReconAction subroutine will calculate the topological charge using the field,
c     strength tensor TopQandReconAction. It is a sum over all elementary 4 plaquettes.
c
c     Modified 26th July 2002 to mean-field improve F_mu_nu by u_0^8
c     Made into a module by F. Bissey in Dec. 2003
c
      MODULE L_TopQandReconAction

      Contains

      subroutine TopQandReconAction(Fr,Fi,uzero,TopQDensity,Q,ReconAction,pluckbar)

      USE GS_LATTICESIZE

      implicit none

      integer,parameter                                       :: nc=3
      integer,parameter                                       :: mu=4
      integer,parameter                                       :: nf=6
      double precision,dimension(nx,ny,nz,nt,nc,nc,nf)        :: Fr,Fi
!HPF$ DISTRIBUTE Fr(*,*,BLOCK,BLOCK,*,*,*)
!HPF$ DISTRIBUTE Fi(*,*,BLOCK,BLOCK,*,*,*)
      double precision                                        :: uzero
      double precision,dimension(nx,ny,nz,nt)                 :: TopQDensity
!HPF$ DISTRIBUTE TopQDensity(*,*,BLOCK,BLOCK)
      double precision                                        :: Q
      double precision,dimension(nx,ny,nz,nt)                 :: ReconAction
!HPF$ DISTRIBUTE ReconAction(*,*,BLOCK,BLOCK)
      double precision                                        :: pluckbar


c     local variables
c
      double precision                                        :: pi
      integer                                                 :: ic,jc,kc,i,j
      integer,dimension(2:mu,1:mu-1)                          :: Findex
!HPF$ DISTRIBUTE Findex(*,*)

      pi = 4.0d0 * atan(1.0d0)

c     Calculate the index that represents the unique combination of
c     mu and nu.
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
c     The unique mu/nu index combination is rolled into the
c     single index by the following relation
c 
      do i=2,mu
         do j =1,i-1

            Findex(i,j) = (i+j) - 3 * ( (i-j) / 3 ) - 1

         end do
      end do

c     TopQDensity(x,y,z,t) is proportional to
c     colourtrace{levicevita(mu,nu,rho,sigma) * F(x,y,z,t,a,b,mu,nu)F(x,y,z,t,a,b,rho,sigma)]}
c

      TopQDensity = 0.0d0

      do ic=1,nc
         do kc=1,nc
 
            TopQDensity(:,:,:,:) = TopQDensity(:,:,:,:) +
     &           ( Fr(:,:,:,:,ic,kc,Findex(2,1)) * Fr(:,:,:,:,kc,ic,Findex(4,3)) -
     &             Fr(:,:,:,:,ic,kc,Findex(3,1)) * Fr(:,:,:,:,kc,ic,Findex(4,2)) +
     &             Fr(:,:,:,:,ic,kc,Findex(4,1)) * Fr(:,:,:,:,kc,ic,Findex(3,2)) -
     &             Fi(:,:,:,:,ic,kc,Findex(2,1)) * Fi(:,:,:,:,kc,ic,Findex(4,3)) +
     &             Fi(:,:,:,:,ic,kc,Findex(3,1)) * Fi(:,:,:,:,kc,ic,Findex(4,2)) -
     &             Fi(:,:,:,:,ic,kc,Findex(4,1)) * Fi(:,:,:,:,kc,ic,Findex(3,2)) ) * 8.0d0

         end do
      end do
c
c  Include mean-field improvement here
c
      TopQDensity(:,:,:,:) = TopQDensity(:,:,:,:) / ( 32.0d0 * ( pi**2 ) * uzero**8 ) 
c
c     Calculation of the topological charge
c
      Q = sum(TopQDensity)
c
c     Now get the reconstructed action from Fmunu
c
      ReconAction = 0.0d0

      do kc = 1,nf
         do ic=1,nc
            do jc=1,nc

               ReconAction(:,:,:,:) = ReconAction(:,:,:,:) +
     &                Fr(:,:,:,:,ic,jc,kc)*Fr(:,:,:,:,jc,ic,kc)
     &              - Fi(:,:,:,:,ic,jc,kc)*Fi(:,:,:,:,jc,ic,kc)

            end do
         end do
      end do

      ReconAction = ReconAction / ( nc * mu*(mu-1) )

      pluckbar = sum(ReconAction) / ( nx*ny*nz*nt )

      return

      end subroutine TopQandReconAction

cccccccccccccccccccccccccccccccccccccccccccccccc

      END MODULE L_TopQandReconAction
