##############################################
#
#  Helix, Orion, APAC and PC settings
#
#COMPILER = pghpf
#CFLAGS   = -I../ -I../../GSnew -Mextend
#helix
#Orion OPT      = -Mmpi -64 -fast -W0,-fast -W0,-xvector=no -xarch=v9a
#Linux OPT      = -byteswapio
#Linux intel fortran
#COMPILER = ifort
#CFLAGS   = -132
#OPT      = -warn all -O2 -ansi_alias -scalar_rep -parallel -static -vec_report5 -convert big_endian
#Linux g95
#COMPILER  = ~/g95-install/bin/i686-unknown-linux-gnu-g95
#CFLAGS    = -ffixed-line-length-132 -ftrace=frame
#OPT       = -Wall -O2 -fendian=big -march=athlon-4
#PREP      = -cpp
#APAC  OPT      = -Mmpi -fast -W0,-fast -Mbyteswapio
#Other OPT      = -byteswapio -fast -W0,-fast -O4 -W0,-O4
COMPILER = gfortran
CFLAGS   = -fconvert=big-endian -ffixed-line-length-132
OPT      = -Wall -O2 -march=native -floop-interchange -floop-strip-mine -floop-block
PREP     = -cpp
LINK     =

.SUFFIXES: .f .o
.f.o:
	$(COMPILER) $(CFLAGS) $(OPT) -c $<

OBJ = gaugetest.o

OBJ2D = gaugetest2d.o

EXECS = GaugeTest \
	GaugeTest2D

all : $(EXECS)

GaugeTest : $(OBJ)
	$(COMPILER) $(CFLAGS) $(OPT) $(PREP) -o GaugeTest $(OBJ) $(LINK)

GaugeTest2D : $(OBJ2D)
	$(COMPILER) $(CFLAGS) $(OPT) $(PREP) -o GaugeTest2D $(OBJ2D) $(LINK)

clean :
	rm -f *.o *~* *.mod $(EXECS)

