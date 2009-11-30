c
c     local variables
c
      integer                                                 		:: it
      integer                                                 		:: ic,jc,kc
c      integer                                                           :: ic1,jc1,kc1
c      integer                                                           :: ic2,jc2,kc2
c
c     diagonal links
c
      double precision,dimension(nx,ny,nz,nt,nc,nc)        		:: rd1r,rd1i,ld1r,ld1i
!HPF$ DISTRIBUTE rd1r(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE rd1i(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE ld1r(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE ld1i(*,*,BLOCK,BLOCK,*,*)
      double precision,dimension(nx,ny,nz,nt,nc,nc)    			:: rd2r,rd2i,ld2r,ld2i
!HPF$ DISTRIBUTE rd2r(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE rd2i(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE ld2r(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE ld2i(*,*,BLOCK,BLOCK,*,*)
      double precision,dimension(nx,ny,nz,nt,nc,nc)        		:: rd3r,rd3i,ld3r,ld3i
!HPF$ DISTRIBUTE rd3r(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE rd3i(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE ld3r(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE ld3i(*,*,BLOCK,BLOCK,*,*)
c
c     bottom links (built up iteratively)
c
      double precision,dimension(nx,ny,nz,nt,nc,nc)           		:: bllr,blli
!HPF$ DISTRIBUTE bllr(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE blli(*,*,BLOCK,BLOCK,*,*)
      double precision,dimension(nx,ny,nz,nt,nc,nc)           		:: brlr,brli
!HPF$ DISTRIBUTE brlr(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE brli(*,*,BLOCK,BLOCK,*,*)
      double precision,dimension(nx,ny,nz,nt,nc,nc)      		:: tlr,tli
!HPF$ DISTRIBUTE tlr(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE tli(*,*,BLOCK,BLOCK,*,*)
      double precision,dimension(nx,ny,nz,nt,nc,nc,nL(that))            :: tltr,tlti
!HPF$ DISTRIBUTE tltr(*,*,BLOCK,BLOCK,*,*,*)
!HPF$ DISTRIBUTE tlti(*,*,BLOCK,BLOCK,*,*,*)
c
c     shifted time links
c
      double precision,dimension(nx,ny,nz,nt,nc,nc)           		:: tl1r,tl1i
!HPF$ DISTRIBUTE tl1r(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE tl1i(*,*,BLOCK,BLOCK,*,*)
      double precision,dimension(nx,ny,nz,nt,nc,nc)           		:: tl2r,tl2i
!HPF$ DISTRIBUTE tl2r(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE tl2i(*,*,BLOCK,BLOCK,*,*)
c
c     top links (time shifted bottom links)
c
      double precision,dimension(nx,ny,nz,nt,nc,nc)                     :: sbllr,sblli
!HPF$ DISTRIBUTE sbllr(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE sblli(*,*,BLOCK,BLOCK,*,*)
      double precision,dimension(nx,ny,nz,nt,nc,nc)                     :: sbrlr,sbrli
!HPF$ DISTRIBUTE sbrlr(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE sbrli(*,*,BLOCK,BLOCK,*,*)
      double precision,dimension(nx,ny,nz,nt,nc,nc)                     :: sllr,slli
!HPF$ DISTRIBUTE sllr(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE slli(*,*,BLOCK,BLOCK,*,*)
      double precision,dimension(nx,ny,nz,nt,nc,nc)                     :: srlr,srli
!HPF$ DISTRIBUTE srlr(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE srli(*,*,BLOCK,BLOCK,*,*)
c
c     temporary arrays
c
      double precision,dimension(nx,ny,nz,nt)                           :: Res
!HPF$ DISTRIBUTE Res(*,*,BLOCK,BLOCK)
      double precision,dimension(nx,ny,nz,nt)                           :: tmp1r,tmp1i
!HPF$ DISTRIBUTE tmp1r(*,*,BLOCK,BLOCK)
!HPF$ DISTRIBUTE tmp1i(*,*,BLOCK,BLOCK)
      double precision,dimension(nx,ny,nz,nt)                           :: tmp2r,tmp2i
!HPF$ DISTRIBUTE tmp2r(*,*,BLOCK,BLOCK)
!HPF$ DISTRIBUTE tmp2i(*,*,BLOCK,BLOCK)
      double precision,dimension(nx,ny,nz,nt)                           :: tmp3r,tmp3i
!HPF$ DISTRIBUTE tmp3r(*,*,BLOCK,BLOCK)
!HPF$ DISTRIBUTE tmp3i(*,*,BLOCK,BLOCK)
      double precision,dimension(nx,ny,nz,nt,nc,nc)                     :: tmpr, tmpi
!HPF$ DISTRIBUTE tmpr(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE tmpi(*,*,BLOCK,BLOCK,*,*)
c
c     epsilon 
c
c      double precision,dimension(nc,nc,nc)                              ::eps
c!HPF$ DISTRIBUTE eps(*,*,*)
