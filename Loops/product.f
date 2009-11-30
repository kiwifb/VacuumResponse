c     subroutine to calculate product of L links in dhat direction
c     made into a module by F. Bissey Dec. 2003
c
      MODULE L_product

      Contains

      subroutine product(ur,ui,dhat,updr,updi,L)

      USE GS_LATTICESIZE

      implicit none

c     global variabes

      integer,parameter                                   :: nc=3 !sigma,colour
      integer,parameter                                   :: mu=4 !directions

      double precision,dimension(nx,ny,nz,nt,mu,nc,nc)    :: ur,ui !links
!HPF$ DISTRIBUTE ur(*,*,BLOCK,BLOCK,*,*,*)
!HPF$ DISTRIBUTE ui(*,*,BLOCK,BLOCK,*,*,*)

      integer                                             :: dhat !direction that the product 
                                                                  !of links are to be calculated in 

      double precision,dimension(nx,ny,nz,nt,nc,nc)       :: updr,updi !product of links in dhat direction
!HPF$ DISTRIBUTE updr(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE updi(*,*,BLOCK,BLOCK,*,*)

      integer                                             :: L !number of links to be included in the product

c     local variables
c     the following temporary arrays hold shifted links and the product of links while the next link is being added on
c
      double precision,dimension(nx,ny,nz,nt,nc,nc)       :: tempr,tempi,usr,usi  
!HPF$ DISTRIBUTE tempr(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE tempi(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE usr  (*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE usi  (*,*,BLOCK,BLOCK,*,*)

      integer                                             :: ic,jc,kc !colour counters

      if (L == 0) then                                              !no links involved in the product, return the identity

                          updr(:,:,:,:,:,:)   = 0.0d0
         forall (ic=1:nc) updr(:,:,:,:,ic,ic) = 1.0d0
                          updi(:,:,:,:,:,:)   = 0.0d0

      else if (L == 1) then !only one link involved in the product

         updr(:,:,:,:,:,:) = ur(:,:,:,:,dhat,:,:)
         updi(:,:,:,:,:,:) = ui(:,:,:,:,dhat,:,:)

      else

c     setting temp equal to the product of L-1 links which were calculated last time product was called
         tempr = updr
         tempi = updi
         usr   = cshift(ur(:,:,:,:,dhat,:,:),dim=dhat,shift=L-1)
         usi   = cshift(ui(:,:,:,:,dhat,:,:),dim=dhat,shift=L-1)
         updr = 0.0d0
         updi = 0.0d0

c     add next link onto temp to give upd
         do ic = 1,nc
            do jc = 1,nc
               do kc = 1,nc

                  updr(:,:,:,:,ic,jc) = updr(:,:,:,:,ic,jc)  
     &                 + tempr(:,:,:,:,ic,kc) * usr(:,:,:,:,kc,jc)
     &                 - tempi(:,:,:,:,ic,kc) * usi(:,:,:,:,kc,jc)
                  updi(:,:,:,:,ic,jc) = updi(:,:,:,:,ic,jc)  
     &                 + tempr(:,:,:,:,ic,kc) * usi(:,:,:,:,kc,jc)
     &                 + tempi(:,:,:,:,ic,kc) * usr(:,:,:,:,kc,jc)

               end do
            end do
         end do

      end if

      return


      end subroutine product

      END MODULE L_product
