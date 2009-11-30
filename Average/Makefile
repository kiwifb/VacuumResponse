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
COMPILER  = ~/g95-install/bin/i686-pc-linux-gnu-g95
CFLAGS    = -I../Loops  -I../GSnew -ffixed-line-length-132
OPT       = -Wall -O2 -march=athlon-4
PREP      = -cpp
#APAC  OPT      = -Mmpi -fast -W0,-fast -Mbyteswapio
#Other OPT      = -byteswapio -fast -W0,-fast -O4 -W0,-O4
LINK     =

.SUFFIXES: .f .o
.f.o:
	$(COMPILER) $(CFLAGS) $(OPT) -c $<

OBJ =   writeAvg.o \
	UniversalAverage.o

BAR   =  ../Loops/baryonParam.o

EXECS = UniversalAverage

all : mbar $(OBJ) $(EXECS)

mbar:
	$(MAKE) -C ../Loops mbar

UniversalAverage : $(OBJ) $(BAR)
	$(COMPILER) $(CFLAGS) $(OPT) $(PREP) -o UniversalAverage $(BAR) $(OBJ) $(LINK)

clean :
	rm -f *.o *~* *.mod $(EXECS)

