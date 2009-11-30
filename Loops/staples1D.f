c     ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     computes the product of links for the action associated 
c       with a link in the that direction.
c
c     author: Frederic Bonnet
c     modified by James Zanotti on 23/3/99 to only calculate staples in 3 spatial dimensions 
c     when smearing in the that direction
c     Made into a module by F. Bissey Dec. 2003
c
      MODULE L_staples1D

      Contains

      subroutine staples1D(ur,ui,stapler,staplei,that,local)

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
      integer                                                 :: that
      logical                                                 :: local
c
c     local variables, temporary product variables
c
      integer,dimension(3)                                    :: yhat
!HPF$ DISTRIBUTE yhat(*)
      double precision,dimension(nx,ny,nz,nt,nc,nc)           :: tr,ti
!HPF$ DISTRIBUTE tr(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE ti(*,*,BLOCK,BLOCK,*,*)
c
      integer                                                 :: ic,jc,kc,ind,i
c
c     starting of the execution commands
c
c     setting up the yhat array
c     the yhat(1)=xhat,yhat(2)=yhat and yhat(3)=zhat when that eq 4
c

      do ind=1,mu-1

         yhat(ind)    = mod( that + ind ,4 )

         if(yhat(ind) == 0) then

            yhat(ind) = 4

         end if

      end do
c
c     calculation of the link products in the positive plaquette.
c
      stapler = 0.0d0
      staplei = 0.0d0
c
c     the tx,ty and tz contour when that eq 4
c
      do i=1,mu-1                     !loop over the yhat array

         tr = 0.0d0
         ti = 0.0d0

         do ic=1,nc
            do jc=1,nc
               do kc=1,nc

                     tr(:,:,:,:,ic,jc) = tr(:,:,:,:,ic,jc)                   +
     &                 ( cshift(ur(:,:,:,:,yhat(i),ic,kc),dim=that,shift=1)  *
     &                   cshift(ur(:,:,:,:,that,jc,kc),dim=yhat(i),shift=1)  +
     &                   cshift(ui(:,:,:,:,yhat(i),ic,kc),dim=that,shift=1)  *
     &                   cshift(ui(:,:,:,:,that,jc,kc),dim=yhat(i),shift=1)  )
                     ti(:,:,:,:,ic,jc) = ti(:,:,:,:,ic,jc)                   +
     &                 ( cshift(ui(:,:,:,:,yhat(i),ic,kc),dim=that,shift=1)  *
     &                   cshift(ur(:,:,:,:,that,jc,kc),dim=yhat(i),shift=1)  -
     &                   cshift(ur(:,:,:,:,yhat(i),ic,kc),dim=that,shift=1)  *
     &                   cshift(ui(:,:,:,:,that,jc,kc),dim=yhat(i),shift=1)  )

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
c     the negative tx,ty and tz contour when that eq 4
c
c     .true. when we calculate the full staple
c     .false. when we calculte plaqbar where
c     only the upper 3 plaquette are needed
c
         if(local) then
c
            tr = 0.0d0
            ti = 0.0d0

            do ic=1,nc
               do jc=1,nc
                  do kc=1,nc

                     tr(:,:,:,:,ic,jc) = tr(:,:,:,:,ic,jc) +
     &                    ( cshift(cshift(ur(:,:,:,:,yhat(i),kc,ic),dim=yhat(i),shift=-1),dim=that,shift=1) *
     &                    cshift(ur(:,:,:,:,that,jc,kc),dim=yhat(i),shift=-1) -
     &                    cshift(cshift(ui(:,:,:,:,yhat(i),kc,ic),dim=yhat(i),shift=-1),dim=that,shift=1)   *
     &                    cshift(ui(:,:,:,:,that,jc,kc),dim=yhat(i),shift=-1) )
                     ti(:,:,:,:,ic,jc) = ti(:,:,:,:,ic,jc) +
     &                    (-cshift(cshift(ur(:,:,:,:,yhat(i),kc,ic),dim=yhat(i),shift=-1),dim=that,shift=1) *
     &                    cshift(ui(:,:,:,:,that,jc,kc),dim=yhat(i),shift=-1) -
     &                    cshift(cshift(ui(:,:,:,:,yhat(i),kc,ic),dim=yhat(i),shift=-1),dim=that,shift=1)   *
     &                    cshift(ur(:,:,:,:,that,jc,kc),dim=yhat(i),shift=-1) )

                  end do
               end do
            end do

            do ic=1,nc
               do jc=1,nc
                  do kc=1,nc

                     stapler(:,:,:,:,ic,jc) = stapler(:,:,:,:,ic,jc) +
     &                    ( tr(:,:,:,:,ic,kc) * cshift(ur(:,:,:,:,yhat(i),kc,jc),dim=yhat(i),shift=-1) -
     &                    ti(:,:,:,:,ic,kc) * cshift(ui(:,:,:,:,yhat(i),kc,jc),dim=yhat(i),shift=-1) )
                     staplei(:,:,:,:,ic,jc) = staplei(:,:,:,:,ic,jc) +
     &                    ( tr(:,:,:,:,ic,kc) * cshift(ui(:,:,:,:,yhat(i),kc,jc),dim=yhat(i),shift=-1) +
     &                    ti(:,:,:,:,ic,kc) * cshift(ur(:,:,:,:,yhat(i),kc,jc),dim=yhat(i),shift=-1) )

                  end do
               end do
            end do
         end if                 !closes the .true. if
      end do                    !closes the i=1,3 loop
    
      return

      end subroutine staples1D

      END MODULE L_staples1D

