c     subroutine to calculate product of L links in dhat direction
c     made into a module by F. Bissey Dec. 2003
c     100117 FB fork to a modules containing more link product functions
c
      MODULE L_product

      Contains
c
c     take lattice links ur and ui a direction dhat, a number of links L
c     and an input/result link variable updr,updi.
c     A product of links in the dhat direction is produced iteratively with
c     function. The link u is shifted L units in the dhat direction then
c     upd is multiplied by the shifted u and the result returned in upd.
c
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
c
c     A more general product function. The previous one can only produce
c     products of links along a line. This function is a 2D generalisation.
c     A product of links in two dimensions is produced iteratively with this function.
c     The link u is shifted ix units in the xdir direction then iy units in
c     the ydir direction. upd is then multiplied by the shifted u in the dhat direction
c     and the result is returned in upd.
c     Note that if dhat is different from xdir and ydir you are going 3D.
c     A 3D version of this routine could be created for such cases.
c
      subroutine product2D(ur,ui,dhat,xdir,ydir,ix,iy,updr,updi)

      USE GS_LATTICESIZE

      implicit none

c     global variabes

      integer,parameter                                   :: nc=3 !sigma,colour
      integer,parameter                                   :: mu=4 !directions

      double precision,dimension(nx,ny,nz,nt,mu,nc,nc)    :: ur,ui !links
!HPF$ DISTRIBUTE ur(*,*,BLOCK,BLOCK,*,*,*)
!HPF$ DISTRIBUTE ui(*,*,BLOCK,BLOCK,*,*,*)

      integer                                             :: dhat,xdir, ydir

      double precision,dimension(nx,ny,nz,nt,nc,nc)       :: updr,updi !product of links
!HPF$ DISTRIBUTE updr(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE updi(*,*,BLOCK,BLOCK,*,*)

      integer                                             :: ix,iy

c     local variables
c     the following temporary arrays hold shifted links and the product of links while the next link is being added on
c
      double precision,dimension(nx,ny,nz,nt,nc,nc)       :: tempr,tempi,usr,usi
!HPF$ DISTRIBUTE tempr(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE tempi(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE usr  (*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE usi  (*,*,BLOCK,BLOCK,*,*)

      integer                                             :: ic,jc,kc !colour counters

      if ((ix == 0) .and. (iy == 0)) then                 !no links involved in the product, return the identity

                          updr(:,:,:,:,:,:)   = 0.0d0
         forall (ic=1:nc) updr(:,:,:,:,ic,ic) = 1.0d0
                          updi(:,:,:,:,:,:)   = 0.0d0

      else if ((ix+iy) == 1) then !only one link involved in the product

         if (ix == 1) then
            updr(:,:,:,:,:,:) = ur(:,:,:,:,xdir,:,:)
            updi(:,:,:,:,:,:) = ui(:,:,:,:,xdir,:,:)
         else
            updr(:,:,:,:,:,:) = ur(:,:,:,:,ydir,:,:)
            updi(:,:,:,:,:,:) = ui(:,:,:,:,ydir,:,:)
         endif
      else

c     setting temp equal to the product of links provide as input
         tempr = updr
         tempi = updi
c     shifting u at (ix,iy) coordinate and setting its direction.
         usr   = cshift(cshift(ur(:,:,:,:,dhat,:,:),dim=xdir,shift=ix),dim=ydir,shift=iy)
         usi   = cshift(cshift(ui(:,:,:,:,dhat,:,:),dim=ydir,shift=ix),dim=ydir,shift=iy)
c     preparing upd for the result
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


      end subroutine products

      END MODULE L_product
