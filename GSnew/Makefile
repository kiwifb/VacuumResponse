##############################################
#
#  Helix, Orion, APAC and PC settings
#
#COMPILER = pghpf
#CFLAGS   = -Mextend
#helix
#OPT	= -Mmpi -O2 -Mnoframe -Mvect=sse -Minfo -byteswapio
#Orion OPT      = -Mmpi -64 -fast -W0,-fast -W0,-xvector=no -xarch=v9a
#Linux OPT      = -byteswapio
#Linux intel fortran
#COMPILER = ifort
#CFLAGS   = -132
#OPT      = -warn all -O2 -ansi_alias -scalar_rep -parallel -static -vec_report5 -convert big_endian
#Linux g95
#COMPILER := ~/g95-v0.9/bin/i686-pc-linux-gnu-g95
COMPILER := ~/g95-install/bin/i686-unknown-linux-gnu-g95
OPT      := -O2 -mtune=athlon-4 -fendian=big -ffixed-line-length-132
UPREP    := -cpp
#APAC  OPT      = -Mmpi -fast -W0,-fast -Mbyteswapio
#Other OPT      = -byteswapio -fast -W0,-fast -O4 -W0,-O4
LINK     :=

.SUFFIXES: .f .o
.f.o:
	$(COMPILER) $(PREP) $(CFLAGS) $(OPT) -c $<

FIRST  := latticeSize.o \
          hpfRandomSeed.o

COMMON := cooling.o \
          fixsu3.o \
          newFixsu3.o \
          parityTransform.o \
          uu.o \
          uu9.o \
          rectangles.o \
          squares.o \
          staples.o \
          getAction.o \
          initialu.o \
          maskWilson.o \
          maskRect.o \
          maskUltraFast.o \
          maskUltra.o \
          setMask.o \
          apeSmear.o \
          pseudoheat.o \
          pseudosweep.o \
          su2random.o \
          su3random.o \
          tadpoleimp.o \
          readLinks.o \
          readLinksRisky.o \
          writeLinks.o

s12: $(FIRST) $(COMMON)
s12: PREP := $(UPREP) -D__S12__

s16: $(FIRST) $(COMMON)
s16: PREP := $(UPREP) -D__S16__

all: s16

clean :
	rm -f *.o *.d *.mod *~* *#* $(EXECS)


##############################################
# CM5 settings
#
# COMPILER = cmf
# CFLAGS   = -cm5 -vu -extend_source -f90syntax
# LINK     = -lcmsslcm5vu
#
# .SUFFIXES: .fcm .o
# .fcm.o:
#        $(COMPILER) $(CFLAGS) -c $<
#
##############################################
