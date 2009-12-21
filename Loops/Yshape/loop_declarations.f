c
c     local variables
c
      integer                                                 :: it
      integer                                                 :: ic,jc,kc
c
c     diagonal links
c
      double precision,dimension(nx,ny,nz,nt,nc,nc)        :: up1x1r,up1x1i,do1x1r,do1x1i
!HPF$ DISTRIBUTE up1x1r(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE up1x1i(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE do1x1r(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE do1x1i(*,*,BLOCK,BLOCK,*,*)
      double precision,dimension(nx,ny,nz,nt,nc,nc)        :: up1x2r,up1x2i,do1x2r,do1x2i
!HPF$ DISTRIBUTE up1x2r(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE up1x2i(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE do1x2r(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE do1x2i(*,*,BLOCK,BLOCK,*,*)
      double precision,dimension(nx,ny,nz,nt,nc,nc)        :: up2x3r,up2x3i,do2x3r,do2x3i
!HPF$ DISTRIBUTE up2x3r(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE up2x3i(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE do2x3r(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE do2x3i(*,*,BLOCK,BLOCK,*,*)
c
c     products of links
c
      double precision,dimension(nx,ny,nz,nt,nc,nc)        :: uptr,upti
!HPF$ DISTRIBUTE uptr(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE upti(*,*,BLOCK,BLOCK,*,*)
c
c     bottom and top links (built up iteratively)
c
      double precision,dimension(nx,ny,nz,nt,nc,nc)        :: bxlr,bxli
!HPF$ DISTRIBUTE bxlr(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE bxli(*,*,BLOCK,BLOCK,*,*)
      double precision,dimension(nx,ny,nz,nt,nc,nc)        :: bllr,blli
!HPF$ DISTRIBUTE bllr(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE blli(*,*,BLOCK,BLOCK,*,*)
      double precision,dimension(nx,ny,nz,nt,nc,nc)        :: brlr,brli
!HPF$ DISTRIBUTE brlr(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE brli(*,*,BLOCK,BLOCK,*,*)
      double precision,dimension(nx,ny,nz,nt,nc,nc)        :: tlr,tli
!HPF$ DISTRIBUTE tlr(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE tli(*,*,BLOCK,BLOCK,*,*)
c
c     shifted time links
c
      double precision,dimension(nx,ny,nz,nt,nc,nc)        :: ustr,usti
!HPF$ DISTRIBUTE ustr(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE usti(*,*,BLOCK,BLOCK,*,*)
c
c     staples/temporary arrays
c
      double precision,dimension(nx,ny,nz,nt,nc,nc)        :: sxlr,sxli
!HPF$ DISTRIBUTE sxlr(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE sxli(*,*,BLOCK,BLOCK,*,*)
      double precision,dimension(nx,ny,nz,nt,nc,nc)        :: sllr,slli
!HPF$ DISTRIBUTE sllr(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE slli(*,*,BLOCK,BLOCK,*,*)
      double precision,dimension(nx,ny,nz,nt,nc,nc)        :: srlr,srli
!HPF$ DISTRIBUTE srlr(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE srli(*,*,BLOCK,BLOCK,*,*)
c
c     temporary arrays
c
      double precision,dimension(nx,ny,nz,nt,nc,nc)        :: tmpr,tmpi
!HPF$ DISTRIBUTE tmpr(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE tmpi(*,*,BLOCK,BLOCK,*,*)
      double precision,dimension(nx,ny,nz,nt)              :: Res
!HPF$ DISTRIBUTE Res(*,*,BLOCK,BLOCK)
