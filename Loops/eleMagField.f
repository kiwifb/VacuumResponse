c     cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     EleMagField subroutine will calculate the sum of the squared 
c        electric and magnetic fields 
c     Made into a module by F. Bissey Dec. 2003
c
      MODULE L_EleMagField

      Contains

      subroutine EleMagField(Fr,Fi,ele,mag)

      USE GS_LATTICESIZE

      implicit none

c     global variables

      integer,parameter                                       :: nc=3
      integer,parameter                                       :: mu=4

      integer,parameter                                       :: nf=6

      double precision,dimension(nx,ny,nz,nt,nc,nc,nf)        :: Fr,Fi
!HPF$ DISTRIBUTE Fr(*,*,BLOCK,BLOCK,*,*,*)
!HPF$ DISTRIBUTE Fi(*,*,BLOCK,BLOCK,*,*,*)
      double precision,dimension(nx,ny,nz,nt)                 :: ele
!HPF$ DISTRIBUTE ele(*,*,BLOCK,BLOCK)
      double precision,dimension(nx,ny,nz,nt)                 :: mag
!HPF$ DISTRIBUTE mag(*,*,BLOCK,BLOCK)

c     local variables

      double precision                                        :: pi
      integer                                                 :: ic,kc,i,j
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

      ele = 0.0d0

      do ic=1,nc
         do kc=1,nc

            ele(:,:,:,:) = ele(:,:,:,:) +
     &           ( Fr(:,:,:,:,ic,kc,Findex(4,3)) * Fr(:,:,:,:,kc,ic,Findex(4,3)) +
     &             Fr(:,:,:,:,ic,kc,Findex(4,2)) * Fr(:,:,:,:,kc,ic,Findex(4,2)) +
     &             Fr(:,:,:,:,ic,kc,Findex(4,1)) * Fr(:,:,:,:,kc,ic,Findex(4,1)) -
     &             Fi(:,:,:,:,ic,kc,Findex(4,3)) * Fi(:,:,:,:,kc,ic,Findex(4,3)) -
     &             Fi(:,:,:,:,ic,kc,Findex(4,2)) * Fi(:,:,:,:,kc,ic,Findex(4,2)) -
     &             Fi(:,:,:,:,ic,kc,Findex(4,1)) * Fi(:,:,:,:,kc,ic,Findex(4,1)) )

         end do
      end do

      mag = 0.0d0

      do ic=1,nc
         do kc=1,nc

            mag(:,:,:,:) = mag(:,:,:,:) +
     &           ( Fr(:,:,:,:,ic,kc,Findex(3,2)) * Fr(:,:,:,:,kc,ic,Findex(3,2)) +
     &             Fr(:,:,:,:,ic,kc,Findex(3,1)) * Fr(:,:,:,:,kc,ic,Findex(3,1)) +
     &             Fr(:,:,:,:,ic,kc,Findex(2,1)) * Fr(:,:,:,:,kc,ic,Findex(2,1)) -
     &             Fi(:,:,:,:,ic,kc,Findex(3,2)) * Fi(:,:,:,:,kc,ic,Findex(3,2)) -
     &             Fi(:,:,:,:,ic,kc,Findex(3,1)) * Fi(:,:,:,:,kc,ic,Findex(3,1)) -
     &             Fi(:,:,:,:,ic,kc,Findex(2,1)) * Fi(:,:,:,:,kc,ic,Findex(2,1)) )

         end do
      end do
c
c  Correct for factor of i in  F_43 = i E_3  etc.
c
      ele = -ele / ( nc )
      mag = mag / ( nc )

      return

      end subroutine EleMagField

      END MODULE L_EleMagField
