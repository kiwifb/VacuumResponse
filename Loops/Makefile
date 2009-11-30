##############################################
#
#  Helix, Orion, APAC and PC settings
#
#COMPILER = pghpf
#CFLAGS   = -Mextend -I../GSnew
#helix
#OPT    = -Mmpi -O2 -Mnoframe -Mvect=sse -Minfo -byteswapio
#LINK   = -L/opt/mpich/gnu/lib
#PREP   = -Mpreprocess
#LINK   = -L/usr/opt/mpi/lib
#Orion OPT      = -Mmpi -64 -fast -W0,-fast -W0,-xvector=no -xarch=v9a
#Linux OPT      = -byteswapio
#Linux intel fortran
#COMPILER = ifort
#CFLAGS   = -132 -I../GSnew
#OPT      = -warn all -O2 -ansi_alias -scalar_rep -parallel -static -vec_report5 -convert big_endian
#Linux g95
COMPILER  := ~/g95-install/bin/i686-unknown-linux-gnu-g95
INCLUDE   := -I../GSnew
CFLAGS    := -ffixed-line-length-132
OPT       := -Wall -O2 -fendian=big -march=athlon-4
UPREP     := -cpp
#APAC  OPT      = -Mmpi -fast -W0,-fast -Mbyteswapio
#Other OPT      = -byteswapio    -fast -W0,-fast   -O4 -W0,-O4
export COMPILER CFLAGS OPT

.SUFFIXES: .f .o
.f.o:
	$(COMPILER) $(INCLUDE) $(CFLAGS) $(OPT) -c $<

COMDIR := ../GSnew

FILES  := VacuumResp.f \
	writeshape.o

FULLCOMMAND := $(COMPILER) $(INCLUDE) $(CFLAGS) $(OPT)

COMMON := $(COMDIR)/fixsu3.o \
          $(COMDIR)/newFixsu3.o \
          $(COMDIR)/readLinks.o \
          $(COMDIR)/getAction.o \
          $(COMDIR)/squares.o \
          $(COMDIR)/staples.o \
          $(COMDIR)/rectangles.o \
          $(COMDIR)/uu.o \
          $(COMDIR)/tadpoleimp.o \
          $(COMDIR)/cooling.o
#         /usr/pgi/src/sunmpi64.o

Dir  :=  epsilonIndex.o \
         product.o \
         staples3D.o \
         ape_smear.o \
         ape_smear3D.o \
         eleMagField.o \
         topQandReconAction.o \
         lotsaLoops.o \
         fMuNu.o \
			boxes.o

BAR   :=  baryonParam.o
#         /usr/pgi/src/sunmpi64.o

Ysh   := Yshape/loops.o

VYsh  := VYshape/loops.o

VY1sh := VY1shape/loops.o

Lsh   := Lshape/loops.o

LPsh   := LPshape/loops.o

TYsh  := TYshape/loops.o

DQsh  := DiQuark/loops.o

QQ    := qqbar/loops.o

DLT   := Delta/loops.o

TEST  := Test/loops.o

EXECS := VacuumRespY \
	VacuumRespL \
	VacuumRespLP \
	VacuumRespTY \
	VacuumRespVY \
	VacuumRespVY1 \
	VacuumRespQQ \
	VacuumRespDQ \
	VacuumRespDLT

TEXECS := VacuumRespTest

all: s16

s16: gsnew16 essential $(EXECS)
s16: PREP := $(UPREP) -D__S16__
	export PREP

s12: gsnew12 essential $(EXECS)
s12: PREP := $(UPREP) -D__S12__
	export PREP

essential : mbar $(Dir)

gsnew16:
	$(MAKE) -C ../GSnew s16

gsnew12:
	$(MAKE) -C ../GSnew s12

mbar: $(BAR)

VacuumRespQQ :$(COMMON) $(BAR) $(Dir)
	$(MAKE) -C qqbar -f ../Templates/Makefile
	$(FULLCOMMAND) -Iqqbar $(PREP) -c writeshape.f
	$(FULLCOMMAND) -Iqqbar $(PREP) -o VacuumRespQQ $(FILES) $(COMMON) $(BAR) $(Dir) $(QQ) $(LINK)

VacuumRespY :$(COMMON) $(BAR) $(Dir)
	$(MAKE) -C Yshape -f ../Templates/Makefile
	$(FULLCOMMAND) -IYshape $(PREP) -c writeshape.f
	$(FULLCOMMAND) -IYshape $(PREP) -o VacuumRespY $(FILES) $(COMMON) $(BAR) $(Dir) $(Ysh) $(LINK)

VacuumRespL :$(COMMON) $(BAR) $(Dir)
	$(MAKE) -C Lshape -f ../Templates/Makefile
	$(FULLCOMMAND) -ILshape $(PREP) -c writeshape.f
	$(FULLCOMMAND) -ILshape $(PREP) -o VacuumRespL $(FILES) $(COMMON) $(BAR) $(Dir) $(Lsh) $(LINK)

VacuumRespLP :$(COMMON) $(BAR) $(Dir)
	$(MAKE) -C LPshape -f ../Templates/Makefile
	$(FULLCOMMAND) -ILPshape $(PREP) -c writeshape.f
	$(FULLCOMMAND) -ILPshape $(PREP) -o VacuumRespLP $(FILES) $(COMMON) $(BAR) $(Dir) $(LPsh) $(LINK)

VacuumRespVY :$(COMMON) $(BAR) $(Dir)
	$(MAKE) -C VYshape -f ../Templates/Makefile
	$(FULLCOMMAND) -IVYshape $(PREP) -c writeshape.f
	$(FULLCOMMAND) -IVYshape $(PREP) -o VacuumRespVY $(FILES) $(COMMON) $(BAR) $(Dir) $(VYsh) $(LINK)

VacuumRespVY1 :$(COMMON) $(BAR) $(Dir)
	$(MAKE) -C VY1shape -f ../Templates/Makefile
	$(FULLCOMMAND) -IVY1shape $(PREP) -c writeshape.f
	$(FULLCOMMAND) -IVY1shape $(PREP) -o VacuumRespVY1 $(FILES) $(COMMON) $(BAR) $(Dir) $(VY1sh) $(LINK)

VacuumRespTY :$(COMMON) $(BAR) $(Dir)
	$(MAKE) -C TYshape -f ../Templates/Makefile
	$(FULLCOMMAND) -ITYshape $(PREP) -c writeshape.f
	$(FULLCOMMAND) -ITYshape $(PREP) -o VacuumRespTY $(FILES) $(COMMON) $(BAR) $(Dir) $(TYsh) $(LINK)

VacuumRespDQ :$(COMMON) $(BAR) $(Dir)
	$(MAKE) -C DiQuark -f ../Templates/Makefile
	$(FULLCOMMAND) -IDiQuark $(PREP) -c writeshape.f
	$(FULLCOMMAND) -IDiQuark $(PREP) -o VacuumRespDQ $(FILES) $(COMMON) $(BAR) $(Dir) $(DQsh) $(LINK)

VacuumRespDLT :$(COMMON) $(BAR) $(Dir)
	$(MAKE) -C Delta -f ../Templates/Makefile
	$(FULLCOMMAND) -IDelta $(PREP) -c writeshape.f
	$(FULLCOMMAND) -IDelta $(PREP) -o VacuumRespDLT $(FILES) $(COMMON) $(BAR) $(Dir) $(DLT) $(LINK)

test : $(COMMON) $(BAR) $(Dir)
	$(MAKE) -C Test -f ../Templates/Makefile
	$(FULLCOMMAND) -ITest $(PREP) -c writeshape.f
	$(FULLCOMMAND) -ITest $(PREP) -o VacuumRespTest $(FILES) $(COMMON) $(BAR) $(Dir) $(TEST) $(LINK)
test : PREP := $(UPREP) -D__S12_

clean :
	rm -f *.o *~* *.mod $(EXECS) $(TEXECS)
	rm -f qqbar/*.o qqbar/*~* qqbar/*.mod
	rm -f Yshape/*.o Yshape/*~* Yshape/*.mod
	rm -f Lshape/*.o Lshape/*~* Lshape/*.mod
	rm -f LPshape/*.o LPshape/*~* LPshape/*.mod
	rm -f VYshape/*.o VYshape/*~* VYshape/*.mod
	rm -f VY1shape/*.o VY1shape/*~* VY1shape/*.mod
	rm -f TYshape/*.o TYshape/*~* TYshape/*.mod
	rm -f DiQuark/*.o DiQuark/*~* DiQuark/*.mod
	rm -f Delta/*.o Delta/*~* Delta/*.mod
	rm -f Test/*.o Test/*~* Test/*.mod
