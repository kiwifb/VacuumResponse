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
COMPILER  = ~/g95-install/bin/i686-unknown-linux-gnu-g95
CFLAGS    = -ffixed-line-length-132 -ftrace=frame
OPT       = -Wall -O2 -fendian=big -march=athlon-4
PREP      = -cpp
#APAC  OPT      = -Mmpi -fast -W0,-fast -Mbyteswapio
#Other OPT      = -byteswapio -fast -W0,-fast -O4 -W0,-O4
LINK     =

.SUFFIXES: .f .o
.f.o:
	$(COMPILER) $(CFLAGS) $(OPT) -c $<

OBJ = gaugetest.o

EXECS = GaugeTest

all : $(OBJ) $(EXECS)

GaugeTest : $(OBJ)
	$(COMPILER) $(CFLAGS) $(OPT) $(PREP) -o GaugeTest $(OBJ) $(LINK)

clean :
	rm -f *.o *~* *.mod $(EXECS)

