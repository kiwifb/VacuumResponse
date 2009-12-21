c
c        creating shifted base
c
      sbxlr = cshift( bxlr,dim=that,shift=it)
      sbxli = cshift( bxli,dim=that,shift=it)

      sbllr = cshift( bllr,dim=that,shift=it)
      sblli = cshift( blli,dim=that,shift=it)

      sbrlr = cshift( brlr,dim=that,shift=it)
      sbrli = cshift( brli,dim=that,shift=it)

      Res = 0.0d0
ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     first index is for the shape of the bottom
c     and the second one is for the shape of the top
c
c     shape description:
c     1)                2)             3)
c           x                x               x
c            \               |\              |
c             \              | \             |
c              x             |  x            |  x
c             /              |               | /
c            /               |               |/
c           x                x               x
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      include'loop11.f'
      include'loop12.f'
      include'loop13.f'
      include'loop21.f'
      include'loop22.f'
      include'loop23.f'
      include'loop31.f'
      include'loop32.f'
      include'loop33.f'
