c
c     local variables
c
      integer                                                 		:: it
      integer                                                 		:: ic,jc,kc
c
c     diagonal links
c
      double precision,dimension(nx,ny,nz,nt,nc,nc)            :: up1x1r,up1x1i,do1x1r,do1x1i
!HPF$ DISTRIBUTE up1x1r(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE up1x1i(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE do1x1r(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE do1x1i(*,*,BLOCK,BLOCK,*,*)
      double precision,dimension(nx,ny,nz,nt,nc,nc)            :: up2x1r,up2x1i,do2x1r,do2x1i
!HPF$ DISTRIBUTE up2x1r(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE up2x1i(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE do2x1r(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE do2x1i(*,*,BLOCK,BLOCK,*,*)
      double precision,dimension(nx,ny,nz,nt,nc,nc)            :: up3x2r,up3x2i,do3x2r,do3x2i
!HPF$ DISTRIBUTE up3x2r(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE up3x2i(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE do3x2r(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE do3x2i(*,*,BLOCK,BLOCK,*,*)
c
c     bottom links (built up iteratively)
c
      double precision,dimension(nx,ny,nz,nt,nc,nc)           		:: bllr,blli
!HPF$ DISTRIBUTE bllr(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE blli(*,*,BLOCK,BLOCK,*,*)
      double precision,dimension(nx,ny,nz,nt,nc,nc)           		:: brlr,brli
!HPF$ DISTRIBUTE brlr(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE brli(*,*,BLOCK,BLOCK,*,*)
      double precision,dimension(nx,ny,nz,nt,nc,nc)               :: bxlr,bxli
!HPF$ DISTRIBUTE bxlr(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE bxli(*,*,BLOCK,BLOCK,*,*)
      double precision,dimension(nx,ny,nz,nt,nc,nc)               :: cxlr,cxli
!HPF$ DISTRIBUTE cxlr(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE cxli(*,*,BLOCK,BLOCK,*,*)
c
c     time links table
c
      double precision,dimension(nx,ny,nz,nt,nc,nc,nL(that))  		:: tlr,tli
!HPF$ DISTRIBUTE tlr(*,*,BLOCK,BLOCK,*,*,*)
!HPF$ DISTRIBUTE tli(*,*,BLOCK,BLOCK,*,*,*)
c
c     top links (time shifted bottom links)
c
      double precision,dimension(nx,ny,nz,nt,nc,nc)                     :: sbllr,sblli
!HPF$ DISTRIBUTE sbllr(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE sblli(*,*,BLOCK,BLOCK,*,*)
      double precision,dimension(nx,ny,nz,nt,nc,nc)                     :: sbrlr,sbrli
!HPF$ DISTRIBUTE sbrlr(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE sbrli(*,*,BLOCK,BLOCK,*,*)
      double precision,dimension(nx,ny,nz,nt,nc,nc)                     :: sbxlr,sbxli
!HPF$ DISTRIBUTE sbxlr(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE sbxli(*,*,BLOCK,BLOCK,*,*)
c
c     staples & kink like links
c     the index 1 to 3 is for the time link involved:
c
c        (1)
c
c              (2)
c
c        (3)
c
      double precision,dimension(nx,ny,nz,nt,nc,nc)                     :: st1r,st1i
!HPF$ DISTRIBUTE st1r(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE st1i(*,*,BLOCK,BLOCK,*,*)
      double precision,dimension(nx,ny,nz,nt,nc,nc)                     :: st2r,st2i
!HPF$ DISTRIBUTE st2r(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE st2i(*,*,BLOCK,BLOCK,*,*)
      double precision,dimension(nx,ny,nz,nt,nc,nc)                     :: st3r,st3i
!HPF$ DISTRIBUTE st3r(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE st3i(*,*,BLOCK,BLOCK,*,*)
c
c     shifted time links - same convention applies
c
      double precision,dimension(nx,ny,nz,nt,nc,nc)               :: t_1r,t_1i
!HPF$ DISTRIBUTE t_1r(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE t_1i(*,*,BLOCK,BLOCK,*,*)
      double precision,dimension(nx,ny,nz,nt,nc,nc)               :: t_2r,t_2i
!HPF$ DISTRIBUTE t_2r(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE t_2i(*,*,BLOCK,BLOCK,*,*)
      double precision,dimension(nx,ny,nz,nt,nc,nc)               :: t_3r,t_3i
!HPF$ DISTRIBUTE t_3r(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE t_3i(*,*,BLOCK,BLOCK,*,*)
c
c     temporary arrays
c
      double precision,dimension(nx,ny,nz,nt)                           :: Res
!HPF$ DISTRIBUTE Res(*,*,BLOCK,BLOCK)
      double precision,dimension(nx,ny,nz,nt,nc,nc)                     :: tmp1r,tmp1i
!HPF$ DISTRIBUTE tmp1r(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE tmp1i(*,*,BLOCK,BLOCK,*,*)
