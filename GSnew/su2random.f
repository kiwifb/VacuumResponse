!
!     cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!     subroutine that calculates the random su2 configuration
!
!     Authors: Frederi! D.R. Bonnet & Derek B. Leinweber
!     Date: June 1998.
!
      MODULE GS_SU2RANDOM

      CONTAINS

      subroutine su2random(ursu2,uisu2)

      USE GS_LATTICESIZE

      implicit none

!     global variables

      integer,parameter                                       :: ncsu2=2
      integer,parameter                                       :: nsigma=ncsu2*ncsu2
      integer,parameter                                       :: mu=4

      double precision,dimension(nx,ny,nz,nt,mu,nsigma)       :: a4vector
!HPF$ DISTRIBUTE a4vector(*,*,BLOCK,BLOCK,*,*)
      double precision,dimension(nx,ny,nz,nt,mu,ncsu2,ncsu2)  :: ursu2,uisu2
!HPF$ DISTRIBUTE ursu2(*,*,BLOCK,BLOCK,*,*,*)
!HPF$ DISTRIBUTE uisu2(*,*,BLOCK,BLOCK,*,*,*)

!     local variables

      integer,parameter                                       :: nr=2
      double precision,parameter                              :: pi=3.141592653d0
      double precision,dimension(nx,ny,nz,nt,mu)              :: norm
!HPF$ DISTRIBUTE norm(*,*,BLOCK,BLOCK,*)
      double precision,dimension(nx,ny,nz,nt,mu,nr)           :: r,theta
!HPF$ DISTRIBUTE r(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE theta(*,*,BLOCK,BLOCK,*,*)

      integer                                                 :: isigma

!     call one set of random numbers and put them in a big array called r,
!     it can be splited into 2, accessing both sides of the array using parameter nr.
!     the random number is in (0,1) taking the nat log maps it onto (-infty,0]
!     multiplying it by -2 maps it to [0,infty)

      CALL RANDOM_NUMBER( r )
      where( r==0.0d0 ) r = 1.0d0
      r = sqrt( -2.0d0*log(r) )
      CALL RANDOM_NUMBER( theta )
      theta = 2.0d0 * pi * theta

!     a4vetor 1 and 2 are independent and gaussian distributed with
!     mean 0 and standard deviation 1. repeat process to get a4ve! 3 and 4.
!     It is the cos and sine that normally distributed.
!     a4vector becomes normaly distributed on S^4

      a4vector(:,:,:,:,:,1) = r(:,:,:,:,:,1) * cos( theta(:,:,:,:,:,1) )
      a4vector(:,:,:,:,:,2) = r(:,:,:,:,:,1) * sin( theta(:,:,:,:,:,1) )
      a4vector(:,:,:,:,:,3) = r(:,:,:,:,:,2) * cos( theta(:,:,:,:,:,2) )
      a4vector(:,:,:,:,:,4) = r(:,:,:,:,:,2) * sin( theta(:,:,:,:,:,2) )

      norm = sqrt( sum( a4vector**2,dim=6 ) )

!     it when we divide by the norm= r1**2+r2**2 that the points are brought
!     back onto the surface of the sphere

      do isigma = 1, nsigma
         a4vector(:,:,:,:,:,isigma) = a4vector(:,:,:,:,:,isigma) / norm(:,:,:,:,:)
      end do

!     converting a4vector to an su2 matrix

      ursu2(:,:,:,:,:,1,1) = a4vector(:,:,:,:,:,4)
      ursu2(:,:,:,:,:,2,2) = a4vector(:,:,:,:,:,4)
      ursu2(:,:,:,:,:,1,2) = a4vector(:,:,:,:,:,2)
      ursu2(:,:,:,:,:,2,1) =-a4vector(:,:,:,:,:,2)
      uisu2(:,:,:,:,:,1,1) = a4vector(:,:,:,:,:,3)
      uisu2(:,:,:,:,:,2,2) =-a4vector(:,:,:,:,:,3)
      uisu2(:,:,:,:,:,1,2) = a4vector(:,:,:,:,:,1)
      uisu2(:,:,:,:,:,2,1) = a4vector(:,:,:,:,:,1)

      return

      end subroutine su2random

      END MODULE GS_SU2RANDOM
