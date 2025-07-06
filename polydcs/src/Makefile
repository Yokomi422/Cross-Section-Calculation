#=================================================================
#                Main *** POLYDCS *** Makefile                   #
#=================================================================
# Please, select your machine architecture in ARCH...
#-----------------------------------------------------------------
#
      ARCH = linux
#
include SYS.$(ARCH)
#
#
# All ARCHs for debugging/profiling thus overriding defaults
#-----------------------------------------------------------------
#   MYFLGS = -O1 -g
#   MYFLGS = -O3 -pg
#
#-----------------------------------------------------------------
# All done! Do not modify below here...
#-----------------------------------------------------------------
    FFLAGS = $(MYFLGS) -I$(INCDIR)
    INCDIR = .
    INSTDIR= ../bin/$(ARCH)
   EXESIZE = '`./calcsize $(INSTDIR)/$@`'
   INSTALL = @echo '===> $@ is built in $(INSTDIR) and requires $(EXESIZE) bytes of RAM'



.SUFFIXES: .o .f .c

.f.o:
	$(FC) $(FFLAGS) -c $*.f
.c.o:
	$(CC) $(CFLAGS) -c $*.c



OBJ =   acborn.o alcalc.o asym.o bkprint.o \
        bornks.o borntel.o borntm.o cross.o gac.o intangr.o \
        lincross.o mat.o mmatrx.o polcross.o read.o tmat.o




default:	nonlinear


nonlinear:	setnl polydcs
setnl:
		/bin/rm -f *.o *.bak core par.h
		/bin/ln -sf INCLUDE/nonlinear.par par.h


linear:		setl polydcsL	
setl:
		/bin/rm -f *.o *.bak core par.h
		/bin/ln -sf INCLUDE/linear.par par.h


clean:
	/bin/rm -f *.o *.bak core par.h

cleanarch:
	/bin/rm -f *.o *.bak core par.h
	/bin/rm -f $(INSTDIR)/polydcsL
	/bin/rm -f $(INSTDIR)/polydcs

cleanall:
	/bin/rm -f *.o *.bak core par.h
	/bin/rm -rf ../bin/*

polydcs:	polydcs.f polydcs.o $(OBJ) par.h
	$(FC) $(FFLAGS) -o $@ polydcs.o $(OBJ) $(LIBS)
	@echo; if [ ! -d ../bin/$(ARCH) ]; then mkdir ../bin/$(ARCH); fi
	@echo; mv polydcs $(INSTDIR)
	$(INSTALL)

polydcsL:	polydcs.f polydcs.o $(OBJ)
	$(FC) $(FFLAGS) -o $@ polydcs.o $(OBJ) $(LIBS)
	@echo; if [ ! -d ../bin/$(ARCH) ]; then mkdir ../bin/$(ARCH); fi
	@echo; mv polydcsL $(INSTDIR)
	$(INSTALL)


