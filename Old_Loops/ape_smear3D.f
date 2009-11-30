c     author: Mark Stanford
c     modified by James Zanotti on 23/3/99 to only allow smearing in the spatial directions
c     converted into a module using other modules on 18/12/03 by F.Bissey
c
      Module L_ape_smear3D

      Contains

      subroutine ape_smear3D(ur,ui,alpha,that)

      USE GS_LATTICESIZE
      USE L_staples3D
      USE GS_Newfixsu3

      implicit none

c  global variables
         
      integer,parameter                                         :: nc=3 !sigma,colour
      integer,parameter                                         :: mu=4 !directions

      double precision,dimension(nx,ny,nz,nt,mu,nc,nc)          :: ur,ui
!HPF$ DISTRIBUTE ur(*,*,BLOCK,BLOCK,*,*,*)
!HPF$ DISTRIBUTE ui(*,*,BLOCK,BLOCK,*,*,*)

      double precision			                        :: alpha ! smearing fraction

c  local variables

      integer                     	                        :: that,ihat
      integer,dimension(3)             	                        :: yhat !array of directions to be smeared

      double precision,dimension(nx,ny,nz,nt,nc,nc)             ::stapler,staplei !smeared links
!HPF$ DISTRIBUTE stapler(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE staplei(*,*,BLOCK,BLOCK,*,*)

      double precision,dimension(nx,ny,nz,nt,mu,nc,nc)          :: ur_prm,ui_prm
!HPF$ DISTRIBUTE ur_prm(*,*,BLOCK,BLOCK,*,*,*)
!HPF$ DISTRIBUTE ui_prm(*,*,BLOCK,BLOCK,*,*,*)

      integer                 		                        :: ic,jc !counters
      double precision			                        :: g ! g = alpha/(2*(mu-2))

c     starting the execution commands
	
      if(that == 1) then
         yhat(1) = 2
         yhat(2) = 3
         yhat(3) = 4
      else if(that == 2) then
         yhat(1) = 3 
         yhat(2) = 1
         yhat(3) = 4
      else if(that == 3) then
         yhat(1) = 1
         yhat(2) = 2
         yhat(3) = 4
      else if(that == 4) then 
         yhat(1) = 1
         yhat(2) = 2
         yhat(3) = 3
      end if

      g = alpha/(2*(mu-2))

c     smear links in the yhat(ihat) direction using only staples in the remaining two yhat directions
      do ihat=1,mu-1

         call staples3D(ur,ui,stapler,staplei,yhat(ihat),that)
         
         do ic=1,nc
            do jc=1,nc
               ur_prm(:,:,:,:,yhat(ihat),ic,jc) = (1-alpha)*ur(:,:,:,:,yhat(ihat),ic,jc) + 
     &              g*stapler(:,:,:,:,jc,ic)
               ui_prm(:,:,:,:,yhat(ihat),ic,jc) = (1-alpha)*ui(:,:,:,:,yhat(ihat),ic,jc) - 
     &              g*staplei(:,:,:,:,jc,ic)

            end do
         end do
      end do 

      ur_prm(:,:,:,:,that,:,:) = ur(:,:,:,:,that,:,:)
      ui_prm(:,:,:,:,that,:,:) = ui(:,:,:,:,that,:,:)

c     return smeared links
c        do ihat = 1,mu-1
c           ur(:,:,:,:,yhat(ihat),:,:) = ur_prm(:,:,:,:,yhat(ihat),:,:)
c           ui(:,:,:,:,yhat(ihat),:,:) = ui_prm(:,:,:,:,yhat(ihat),:,:)
c        end do


      call Newfixsu3(ur,ui,ur_prm,ui_prm,8)

      end subroutine ape_smear3D 

      End Module L_ape_smear3D
