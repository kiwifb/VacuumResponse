        integer,parameter                                       :: nc=3         !sigma,colour
        integer,parameter                                       :: mu=4         !directions
        integer,parameter                                       :: xhat=4       !spatial direction for wilson loop
        integer,parameter                                       :: yhat=3       !spatial direction for wilson loop
        integer,parameter                                       :: zhat=2       !spatial direction for wilson loop
        integer,parameter                                       :: that=1       !time direction for wilson loop
c        integer,dimension(mu)                                   :: nL = (/ 8, 5*ny/8, 5*nz/8, 5*nt/8 /)
        integer,parameter,dimension(mu)                         :: nL = (/ 8, ny/2 - 1, nz/2 - 1, nt/2 - 2 /)

