##############################################
#
#  Helix, Orion, APAC and PC settings
#
#COMPILER = pghpf
#CFLAGS   = -Mextend -I../GSnew
#helix
#OPT    = -Mmpi -O2 -Mnoframe -Mvect=sse -Minfo -byteswapio 
#LINK   = -L/opt/mpich/gnu/lib
#LINK   = -L/usr/opt/mpi/lib 
#Orion OPT      = -Mmpi -64 -fast -W0,-fast -W0,-xvector=no -xarch=v9a 
#Linux OPT      = -byteswapio 
#Linux intel fortran
#COMPILER = ifort
#CFLAGS   = -132 -I../GSnew
#OPT      = -warn all -O2 -ansi_alias -scalar_rep -parallel -static -vec_report5 -convert big_endian
#Linux g95
COMPILER  = ~/g95-v0.9/bin/i686-pc-linux-gnu-g95
CFLAGS    = -I../GSnew
OPT       = -O2 -mtune=athlon-4 -fendian=big -ffixed-line-length-132
#APAC  OPT      = -Mmpi -fast -W0,-fast -Mbyteswapio
#Other OPT      = -byteswapio    -fast -W0,-fast   -O4 -W0,-O4

.SUFFIXES: .f .o 
.f.o:
	$(COMPILER) $(CFLAGS) $(OPT) -c $< 

COMDIR = ../GSnew

COMMON = $(COMDIR)/fixsu3.o \
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

Dir   =  epsilonIndex.o \
         product.o \
         staples3D.o \
         ape_smear.o \
         ape_smear3D.o \
         eleMagField.o \
         topQandReconAction.o \
         lotsaLoops.o \
         fMuNu.o 

BAR   =  baryonParam.o 
#         /usr/pgi/src/sunmpi64.o

Avg   = AverageLT.o

Yall  = writeAvg.o \
        AverageYall.o 

Tall  = writeAvg.o \
        AverageTall.o

#lop   = lt_loops.o \
#	writeLTshape.o \
#        VacuumRespLT_loop.o

#plan  = lt_loops.o \
#	writeLTshape.o \
#        VacuumRespLT_plan.o

lopY  = y_loops.o \
	writeYshape.o \
	VacuumRespY_loop.o

planY = y_loops.o \
	writeYshape.o \
        VacuumRespY_plan.o

#planD = delta_loops.o \
#	writeDshape.o \
#        VacuumRespD_plan.o

Lplan = ly_loop_rx.o \
	writeLYshape.o \
	VacuumRespLY_plan.o 

Tplan = t_loops.o \
	writeLTshape.o \
	VacuumRespT_plan.o

TYplan= ty_loops.o \
	writeYshape.o \
	VacuumRespTY_plan.o

DQplan= dq_loops.o \
	writeDQshape.o \
	VacuumRespTdq_plan.o

QQ    = qq_loops.o \
	writeDQshape.o \
	writeAvg.o \
	VacuumRespQQ_plan.o

#plY89 = y_loops89.o \
#	writeYshape89.o \
#	VacuumRespY-plan89.o
#
#merge = merge89.o

EXECS = AverageLT \
        AverageYall \
	AverageTall \
        VacuumRespYloop \
        VacuumRespYplan \
	VacuumRespTplan \
	VacuumRespTYplan \
	VacuumRespLYplan \
	VacuumRespDQplan \
	VacuumRespQQ
#	VacuumRespDplan \

all: $(OBJS) $(EXECS)


AverageLT : $(COMMON) $(BAR) $(Dir) $(Avg)
	$(COMPILER) $(CFLAGS) $(OPT) -o AverageLT $(COMMON) $(BAR) $(Avg) $(LINK)

AverageYall : $(COMMON) $(BAR) $(Dir) $(Yall)
	$(COMPILER) $(CFLAGS) $(OPT) -o AverageYall $(COMMON) $(BAR) $(Yall) $(LINK)

AverageTall : $(COMMON) $(BAR) $(Dir) $(Tall)
	$(COMPILER) $(CFLAGS) $(OPT) -o AverageTall $(COMMON) $(BAR) $(Tall) $(LINK)

#VacuumRespLTloop :$(COMMON) $(BAR) $(Dir) $(lop)
#	$(COMPILER) $(CFLAGS) $(OPT) -o VacuumRespLTloop $(COMMON) $(BAR) $(Dir) $(lop) $(LINK)

#VacuumRespLTplan :$(COMMON) $(BAR) $(Dir) $(plan)
#	$(COMPILER) $(CFLAGS) $(OPT) -o VacuumRespLTplan $(COMMON) $(BAR) $(Dir) $(plan) $(LINK)

VacuumRespYloop :$(COMMON) $(BAR) $(Dir) $(lopY)
	$(COMPILER) $(CFLAGS) $(OPT) -o VacuumRespYloop $(COMMON) $(BAR) $(Dir) $(lopY) $(LINK)	

VacuumRespYplan :$(COMMON) $(BAR) $(Dir) $(planY)
	$(COMPILER) $(CFLAGS) $(OPT) -o VacuumRespYplan $(COMMON) $(BAR) $(Dir) $(planY) $(LINK)

#VacuumRespDplan :$(COMMON) $(BAR) $(Dir) $(planD)
#	$(COMPILER) $(CFLAGS) $(OPT) -o VacuumRespDplan $(COMMON) $(BAR) $(Dir) $(planD) $(LINK)

VacuumRespLYplan :$(COMMON) $(BAR) $(Dir) $(Lplan)
	$(COMPILER) $(CFLAGS) $(OPT) -o VacuumRespLYplan $(COMMON) $(BAR) $(Dir) $(Lplan) $(LINK)

VacuumRespTplan :$(COMMON) $(BAR) $(Dir) $(Tplan)
	$(COMPILER) $(CFLAGS) $(OPT) -o VacuumRespTplan $(COMMON) $(BAR) $(Dir) $(Tplan) $(LINK)

VacuumRespTYplan :$(COMMON) $(BAR) $(Dir) $(TYplan)
	$(COMPILER) $(CFLAGS) $(OPT) -o VacuumRespTYplan $(COMMON) $(BAR) $(Dir) $(TYplan) $(LINK)

VacuumRespDQplan :$(COMMON) $(BAR) $(Dir) $(DQplan)
	$(COMPILER) $(CFLAGS) $(OPT) -o VacuumRespDQplan $(COMMON) $(BAR) $(Dir) $(DQplan) $(LINK)

VacuumRespQQ :$(COMMON) $(BAR) $(Dir) $(QQ)
	$(COMPILER) $(CFLAGS) $(OPT) -o VacuumRespQQ $(COMMON) $(BAR) $(Dir) $(QQ) $(LINK)

#VacuumRespYplan89 :$(COMMON) $(Dir) $(plY89)
#	$(COMPILER) $(CFLAGS) $(OPT) -o VacuumRespYplan89 $(COMMON) $(Dir) $(plY89) $(LINK)

#merge89 :$(COMMON) $(Dir) $(merge)
#	$(COMPILER) $(CFLAGS) $(OPT) -o merge89 $(COMMON) $(Dir) $(merge) $(LINK)

clean :
	rm -f *.o *~* *.mod *#* $(EXECS)
