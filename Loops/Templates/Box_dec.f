!
!  universal boxes declarations and initial zeroing.
!
      double precision,dimension(nx,ny,nz,nt,mu,nc,nc)        :: ur,ui !links
!HPF$ DISTRIBUTE ur(*,*,BLOCK,BLOCK,*,*,*)
!HPF$ DISTRIBUTE ui(*,*,BLOCK,BLOCK,*,*,*)
      double precision,dimension(nx,ny,nz,nt,nc,nc)           :: boxr,boxi
!HPF$ DISTRIBUTE boxr(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE boxi(*,*,BLOCK,BLOCK,*,*)

      double precision,dimension(nx,ny,nz,nt,nc,nc)           :: ux1r,ux1i,ux2r,ux2i,uy1r,uy1i,uy2r,uy2i
!HPF$ DISTRIBUTE ux1r(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE ux1i(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE ux2r(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE ux2i(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE uy1r(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE uy1i(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE uy2r(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE uy2i(*,*,BLOCK,BLOCK,*,*)

      integer                                                 :: xdir,ydir
      integer                                                 :: ic,jc,kc

      boxr = 0.d0
      boxi = 0.d0
