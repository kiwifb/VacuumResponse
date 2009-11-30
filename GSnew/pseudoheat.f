!
!     cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!     subroutine that implements the pseudo-heatbath algorithm
!
!     Author: Frederi! D.R. Bonnet & Derek B. Leinweber
!     Date: 28th of July 1998
!
      MODULE GS_PSEUDOHEAT

      CONTAINS

      subroutine pseudoheat(urnewsu2,uinewsu2,phbsr,phbsi,mask,imask,beta,ihat)

      USE GS_LATTICESIZE

      implicit none

!     global variables

      integer,parameter                                       :: ncsu2=2
      integer,parameter                                       :: nsigma=ncsu2*ncsu2
      integer,parameter                                       :: mu=4
      integer,parameter                                       :: nmask=16

      integer                                                 :: ihat,imask
      double precision                                        :: beta

      double precision,dimension(nx,ny,nz,nt,ncsu2,ncsu2)     :: urnewsu2,uinewsu2
      double precision,dimension(nx,ny,nz,nt,ncsu2,ncsu2)     :: phbsr,phbsi
!HPF$ DISTRIBUTE urnewsu2(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE uinewsu2(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE phbsr(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE phbsi(*,*,BLOCK,BLOCK,*,*)
      logical,dimension(nx,ny,nz,nt,mu,nmask)                 :: mask
!HPF$ DISTRIBUTE mask(*,*,BLOCK,BLOCK,*,*)

!     local variables

      integer,parameter                                       :: nr=2
      double precision,parameter                              :: pi=3.141592653d0

      double precision,dimension(nx,ny,nz,nt)                 :: k,x,probrjk,random
!HPF$ DISTRIBUTE k(*,*,BLOCK,BLOCK)
!HPF$ DISTRIBUTE x(*,*,BLOCK,BLOCK)
!HPF$ DISTRIBUTE probrjk(*,*,BLOCK,BLOCK)
!HPF$ DISTRIBUTE random(*,*,BLOCK,BLOCK)
      double precision,dimension(nx,ny,nz,nt)                 :: norm
!HPF$ DISTRIBUTE norm(*,*,BLOCK,BLOCK)
      logical,dimension(nx,ny,nz,nt)                          :: update
!HPF$ DISTRIBUTE update(*,*,BLOCK,BLOCK)
      double precision,dimension(nx,ny,nz,nt,nr)              :: r,theta
!HPF$ DISTRIBUTE r(*,*,BLOCK,BLOCK,*)
!HPF$ DISTRIBUTE theta(*,*,BLOCK,BLOCK,*)
      double precision,dimension(nx,ny,nz,nt,nsigma)          :: a4vector
!HPF$ DISTRIBUTE a4vector(*,*,BLOCK,BLOCK,*)
      double precision,dimension(nx,ny,nz,nt,ncsu2,ncsu2)     :: urprmsu2,uiprmsu2
!HPF$ DISTRIBUTE urprmsu2(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE uiprmsu2(*,*,BLOCK,BLOCK,*,*)

      integer                                                 :: isigma,ic,jc,kc
      integer                                                 :: counter

!     start of the execution commands

      CALL RANDOM_NUMBER( r )
      where( r==0.0d0 ) r = 1.0d0
      r = sqrt( -2.0d0 * log(r) )
      CALL RANDOM_NUMBER( theta )
      theta = 2.0d0 * pi * theta

      a4vector(:,:,:,:,1) = r(:,:,:,:,1) * cos( theta(:,:,:,:,1) )
      a4vector(:,:,:,:,2) = r(:,:,:,:,1) * sin( theta(:,:,:,:,1) )
      a4vector(:,:,:,:,3) = r(:,:,:,:,2) * cos( theta(:,:,:,:,2) )
      a4vector(:,:,:,:,4) = 0.0d0

!     calculates sqrt(a1^2+a2^2+a3^2+a4^2)=norm

      norm = sqrt( sum( a4vector**2,dim=5 ) )

!     calculates the determinant k=|\sum_{\alpha=1}^6\widetilde{U}_{\alpha}|^{1/2}
!     where $\widetilde{U}_{\alpha}\equiv$ the six product of the three links
!     variable which interact with the link in question, i.e. stapler and staplei

        k = sqrt( abs( phbsr(:,:,:,:,1,1) * phbsr(:,:,:,:,2,2) -
     &                 phbsr(:,:,:,:,1,2) * phbsr(:,:,:,:,2,1) -
     &                 phbsi(:,:,:,:,1,1) * phbsi(:,:,:,:,2,2) +
     &                 phbsi(:,:,:,:,1,2) * phbsi(:,:,:,:,2,1) ) )

!     the values of update is all true when imask=1 and all false when imask=2
!     the value of imask is passed in in the subroutine list and defined in the
!     subroutine pseudosweep.

      update = mask(:,:,:,:,ihat,imask)
      counter = 0

!     the do while makes sure that every links are updated on every sweep.
!     an array dimensioned according to the volume of the lattice is reshaped
!     and store in the variable x. x is a random double precision number, uniformely distributed
!     within the region exp(-2*beta*k)<x<1.
!     This generates a4vector distributed with exponential weight
!     exp(beta*k*a4vector).
!     To summarize, the algorithm begins with a trial a4vector=1+ln(x)/(beta*k)
!     where x is a random number uniformely distributed in the region
!     exp(-2*beta*k)<x<1. This generates a4vector distributed with exponential weight
!     exp(beta*k*a4vector). To correct for the factor (1-a4vector**2)^1/2 in
!     P(a4vector)~(1-a4vector**2)^1/2*exp(beta*k*a4vector), reject this a4vector
!     with probability 1.0d0 - sqrt( 1.0d0 - a4vector**2 ) and select a new trial
!     a4vector, repeat this until an a4vector is accepted.

!     update starts to be true if probrjk is false then we have a true with a
!     false which is false then update = ( update .and. probrjk >= random ) becomes
!     true with a true which is true then we go around the while loop again. This
!     time update is still true but probrjk < random is true then the where
!     statement is true and a4vector=x (the trial). This means that
!     update = ( update .and. probrjk >= random )
!     becomes true with false which is false, therefore we exit the while loop.

      do while( any(update) )
         counter = counter + 1
         call RANDOM_NUMBER( x )
         x = ( 1.0d0 - exp( -2.0d0 * beta * k ) ) * x + exp( -2.0d0 * beta * k )
         x = 1.0d0 + log( x ) / ( beta * k )
         probrjk = 1.0d0 - sqrt( 1.0d0 - x**2 )
         call RANDOM_NUMBER( random )
         where( update .and. probrjk < random ) a4vector(:,:,:,:,4) = x(:,:,:,:)
         update = ( update .and. probrjk >= random )
      end do

      do isigma=1,nsigma-1
        a4vector(:,:,:,:,isigma) =
     &     ( a4vector(:,:,:,:,isigma) * sqrt( 1.0d0 - a4vector(:,:,:,:,4)**2 ) )
     &     / norm(:,:,:,:)
      end do

!     this calculates the link variable U=a4vector(4)*I+ia4vector(1to3)*sigma
!     as in su2random, except here we are only interested in one specific
!     direction ihat, the value of ihat is defined according to the loop over
!     the 4 directions imu=1,mu in the subroutine pseudosweep.

!     converting a4vector to an su2 matrix

      urprmsu2(:,:,:,:,1,1) = a4vector(:,:,:,:,4)
      urprmsu2(:,:,:,:,2,2) = a4vector(:,:,:,:,4)
      urprmsu2(:,:,:,:,1,2) = a4vector(:,:,:,:,2)
      urprmsu2(:,:,:,:,2,1) =-a4vector(:,:,:,:,2)
      uiprmsu2(:,:,:,:,1,1) = a4vector(:,:,:,:,3)
      uiprmsu2(:,:,:,:,2,2) =-a4vector(:,:,:,:,3)
      uiprmsu2(:,:,:,:,1,2) = a4vector(:,:,:,:,1)
      uiprmsu2(:,:,:,:,2,1) = a4vector(:,:,:,:,1)

!     this calculates U --> U'=U * staple^{\dag} / k, the U (full link) coming out
!     of su2random is here replaced by U'(partial link, depends on ihat the direction)

      urnewsu2 = 0.0d0
      uinewsu2 = 0.0d0
      do ic=1,ncsu2
         do jc=1,ncsu2
            do kc=1,ncsu2
                  urnewsu2(:,:,:,:,ic,jc) = urnewsu2(:,:,:,:,ic,jc) +
     &              ( urprmsu2(:,:,:,:,ic,kc) * phbsr(:,:,:,:,jc,kc) +
     &                uiprmsu2(:,:,:,:,ic,kc) * phbsi(:,:,:,:,jc,kc) )
     &              / k(:,:,:,:)
                  uinewsu2(:,:,:,:,ic,jc) = uinewsu2(:,:,:,:,ic,jc) +
     &              ( uiprmsu2(:,:,:,:,ic,kc) * phbsr(:,:,:,:,jc,kc) -
     &                urprmsu2(:,:,:,:,ic,kc) * phbsi(:,:,:,:,jc,kc) )
     &              / k(:,:,:,:)
            end do
         end do
      end do

      return

      end subroutine pseudoheat

      END MODULE GS_PSEUDOHEAT
