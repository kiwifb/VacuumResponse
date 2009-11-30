c    author: Frederic Bonnet
c     modified by James Zanotti on 23/3/99 to only calculate staples in 2 spatial dimensions 
c     when smearing in the other spatial direction
c     Made into a module by F. Bissey Dec. 2003
c
      MODULE L_staples3D

      Contains

      subroutine staples3D(ur,ui,stapler,staplei,xhat,that)

      USE GS_LATTICESIZE

      implicit none
c
c     global variable
c
      integer,parameter                                       :: nc=3
      integer,parameter                                       :: mu=4
c
      double precision,dimension(nx,ny,nz,nt,nc,nc)           :: stapler,staplei
!HPF$ DISTRIBUTE stapler(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE staplei(*,*,BLOCK,BLOCK,*,*)
      double precision,dimension(nx,ny,nz,nt,mu,nc,nc)        :: ur,ui
!HPF$ DISTRIBUTE ur(*,*,BLOCK,BLOCK,*,*,*)
!HPF$ DISTRIBUTE ui(*,*,BLOCK,BLOCK,*,*,*)
c
      integer                                                 :: xhat,that,ihat
c
c     local variables, temporary product variables
c
      integer,dimension(2)                                    :: yhat
!HPF$ DISTRIBUTE yhat(*)
      double precision,dimension(nx,ny,nz,nt,nc,nc)           :: tr,ti
!HPF$ DISTRIBUTE tr(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE ti(*,*,BLOCK,BLOCK,*,*)
c
      integer                                                 :: ic,jc,kc,i
      integer                                                 :: iy
c
c     this subroutine differs from the original staples only through the 
c     next block of code; it reads in the value of that; and the staples 
c     are only calculated in 1 or 2 directions rather than 3.
c
c     starting of the execution commands
c
c     setting up the yhat array
c

      iy = 0

      do ihat = 1,4

         if(ihat /= xhat .and. ihat /= that) then

            iy = iy +1
            yhat(iy) = ihat

         end if

      end do
c
c     calculation of the link products in the positive plaquette.
c
      stapler = 0.0d0
      staplei = 0.0d0
c
c
      do i=1,mu-2                     !loop over the yhat array

         tr = 0.0d0
         ti = 0.0d0

         do ic=1,nc
            do jc=1,nc
               do kc=1,nc

                     tr(:,:,:,:,ic,jc) = tr(:,:,:,:,ic,jc)                   +
     &                 ( cshift(ur(:,:,:,:,yhat(i),ic,kc),dim=xhat,shift=1)  *
     &                   cshift(ur(:,:,:,:,xhat,jc,kc),dim=yhat(i),shift=1)  +
     &                   cshift(ui(:,:,:,:,yhat(i),ic,kc),dim=xhat,shift=1)  *
     &                   cshift(ui(:,:,:,:,xhat,jc,kc),dim=yhat(i),shift=1)  )
                     ti(:,:,:,:,ic,jc) = ti(:,:,:,:,ic,jc)                   +
     &                 ( cshift(ui(:,:,:,:,yhat(i),ic,kc),dim=xhat,shift=1)  *
     &                   cshift(ur(:,:,:,:,xhat,jc,kc),dim=yhat(i),shift=1)  -
     &                   cshift(ur(:,:,:,:,yhat(i),ic,kc),dim=xhat,shift=1)  *
     &                   cshift(ui(:,:,:,:,xhat,jc,kc),dim=yhat(i),shift=1)  )

               end do
            end do
         end do

         do ic=1,nc
            do jc=1,nc
               do kc=1,nc

                     stapler(:,:,:,:,ic,jc) = stapler(:,:,:,:,ic,jc)    +
     &                 ( tr(:,:,:,:,ic,kc) * ur(:,:,:,:,yhat(i),jc,kc)  +
     &                   ti(:,:,:,:,ic,kc) * ui(:,:,:,:,yhat(i),jc,kc)  )
                     staplei(:,:,:,:,ic,jc) = staplei(:,:,:,:,ic,jc)   +
     &                 ( ti(:,:,:,:,ic,kc) * ur(:,:,:,:,yhat(i),jc,kc)  -
     &                   tr(:,:,:,:,ic,kc) * ui(:,:,:,:,yhat(i),jc,kc)  )

               end do
            end do
         end do
c
c     calculation of the link products in the negative plaquette
c     starting point in the x direction
c
         tr = 0.0d0
         ti = 0.0d0

         do ic=1,nc
            do jc=1,nc
               do kc=1,nc

                     tr(:,:,:,:,ic,jc) = tr(:,:,:,:,ic,jc)                         +
     &       ( cshift(cshift(ur(:,:,:,:,yhat(i),kc,ic),dim=yhat(i),shift=-1),dim=xhat,shift=1) *
     &                   cshift(ur(:,:,:,:,xhat,jc,kc),dim=yhat(i),shift=-1)          -
     &         cshift(cshift(ui(:,:,:,:,yhat(i),kc,ic),dim=yhat(i),shift=-1),dim=xhat,shift=1) *
     &                   cshift(ui(:,:,:,:,xhat,jc,kc),dim=yhat(i),shift=-1) )
                     ti(:,:,:,:,ic,jc) = ti(:,:,:,:,ic,jc)                         +
     &       (-cshift(cshift(ur(:,:,:,:,yhat(i),kc,ic),dim=yhat(i),shift=-1),dim=xhat,shift=1) *
     &                   cshift(ui(:,:,:,:,xhat,jc,kc),dim=yhat(i),shift=-1)          -
     &         cshift(cshift(ui(:,:,:,:,yhat(i),kc,ic),dim=yhat(i),shift=-1),dim=xhat,shift=1) *
     &                   cshift(ur(:,:,:,:,xhat,jc,kc),dim=yhat(i),shift=-1) )

               end do
            end do
         end do

         do ic=1,nc
            do jc=1,nc
               do kc=1,nc
                     stapler(:,:,:,:,ic,jc) = stapler(:,:,:,:,ic,jc)            +
     &       ( tr(:,:,:,:,ic,kc) * cshift(ur(:,:,:,:,yhat(i),kc,jc),dim=yhat(i),shift=-1) -
     &         ti(:,:,:,:,ic,kc) * cshift(ui(:,:,:,:,yhat(i),kc,jc),dim=yhat(i),shift=-1) )
                     staplei(:,:,:,:,ic,jc) = staplei(:,:,:,:,ic,jc)            +
     &       ( tr(:,:,:,:,ic,kc) * cshift(ui(:,:,:,:,yhat(i),kc,jc),dim=yhat(i),shift=-1) +
     &         ti(:,:,:,:,ic,kc) * cshift(ur(:,:,:,:,yhat(i),kc,jc),dim=yhat(i),shift=-1) )
               end do
            end do
         end do
      end do                                                !closes the i=1,ie loop

      return

      end subroutine staples3D

      END MODULE L_staples3D
