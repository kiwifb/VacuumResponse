c
c     local variables
c
      integer                                                 		:: it
      integer                                                 		:: ic,jc,kc
c
c     diagonal links
c
      double precision,dimension(nx,ny,nz,nt,nc,nc)            :: up1x1r,up1x1i
!HPF$ DISTRIBUTE up1x1r(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE up1x1i(*,*,BLOCK,BLOCK,*,*)
      double precision,dimension(nx,ny,nz,nt,nc,nc)            :: up2x1r,up2x1i
!HPF$ DISTRIBUTE up2x1r(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE up2x1i(*,*,BLOCK,BLOCK,*,*)
      double precision,dimension(nx,ny,nz,nt,nc,nc)            :: up3x2r,up3x2i
!HPF$ DISTRIBUTE up3x2r(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE up3x2i(*,*,BLOCK,BLOCK,*,*)
c
c     bottom links (built up iteratively)
c
      double precision,dimension(nx,ny,nz,nt,nc,nc)           		:: bllr,blli
!HPF$ DISTRIBUTE bllr(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE blli(*,*,BLOCK,BLOCK,*,*)
      double precision,dimension(nx,ny,nz,nt,nc,nc)           		:: bxlr,bxli
!HPF$ DISTRIBUTE brlr(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE brli(*,*,BLOCK,BLOCK,*,*)
      double precision,dimension(nx,ny,nz,nt,nc,nc)               :: cxlr,cxli
!HPF$ DISTRIBUTE cxlr(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE cxli(*,*,BLOCK,BLOCK,*,*)
      double precision,dimension(nx,ny,nz,nt,nc,nc,nL(that))  		:: tlr,tli
!HPF$ DISTRIBUTE tlr(*,*,BLOCK,BLOCK,*,*,*)
!HPF$ DISTRIBUTE tli(*,*,BLOCK,BLOCK,*,*,*)
c
c     shifted time links
c
      double precision,dimension(nx,ny,nz,nt,nc,nc)           		:: t_tr,t_ti
!HPF$ DISTRIBUTE t_tr(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE t_ti(*,*,BLOCK,BLOCK,*,*)
      double precision,dimension(nx,ny,nz,nt,nc,nc)           		:: t_ur,t_ui
!HPF$ DISTRIBUTE t_ur(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE t_ui(*,*,BLOCK,BLOCK,*,*)
      double precision,dimension(nx,ny,nz,nt,nc,nc)               :: t_dr,t_di
!HPF$ DISTRIBUTE t_dr(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE t_di(*,*,BLOCK,BLOCK,*,*)
c
c     top links (time shifted bottom links)
c
      double precision,dimension(nx,ny,nz,nt,nc,nc)                     :: sbllr,sblli
!HPF$ DISTRIBUTE sbllr(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE sblli(*,*,BLOCK,BLOCK,*,*)
      double precision,dimension(nx,ny,nz,nt,nc,nc)                     :: sbxlr,sbxli
!HPF$ DISTRIBUTE sbrlr(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE sbrli(*,*,BLOCK,BLOCK,*,*)
c
c     staples
c
      double precision,dimension(nx,ny,nz,nt,nc,nc)                     :: sllr,slli
!HPF$ DISTRIBUTE sbllr(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE sblli(*,*,BLOCK,BLOCK,*,*)
      double precision,dimension(nx,ny,nz,nt,nc,nc)                     :: sxlr,sxli
!HPF$ DISTRIBUTE sbrlr(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE sbrli(*,*,BLOCK,BLOCK,*,*)
c
c     temporary arrays
c
      double precision,dimension(nx,ny,nz,nt)                           :: Res
!HPF$ DISTRIBUTE Res(*,*,BLOCK,BLOCK)
      double precision,dimension(nx,ny,nz,nt,nc,nc)                     :: tmp1r,tmp1i
!HPF$ DISTRIBUTE tmp1r(*,*,BLOCK,BLOCK,*,*)
!HPF$ DISTRIBUTE tmp1i(*,*,BLOCK,BLOCK,*,*)
