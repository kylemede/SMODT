.c.o:
	$(CC) -c $(CFLAGS) $*.c

CC = g++ -std=c++11
CFLAGS = -Wall -pedantic -fpermissive -g -Werror
LDFLAGS = -pg
LIBS = -lm -lrt 

#---------------------------------------------------------------
PROGSthi = THItester
OBJSthi = Toolboxes/DItools.o THItester.o Toolboxes/OrbSimGeneralToolbox.o

THItester: $(OBJSthi)
	$(CC) -o THItester $(OBJSthi) $(LDFLAGS) $(LIBS)

#---------------------------------------------------------------
PROGSmc = mcONLYorbSimulator
OBJSmc = mcONLYorbSimulator.o rnd/sfmt.o Toolboxes/orbSimStructures.o Toolboxes/RVtools.o Toolboxes/DItools.o Toolboxes/OrbSimGeneralToolbox.o Toolboxes/SimSettingsObj.o Toolboxes/DataObj.o

mcONLYorbSimulator: $(OBJSmc)
	$(CC) -o mcONLYorbSimulator $(OBJSmc) $(LDFLAGS) $(LIBS)
	
#---------------------------------------------------------------

PROGSsa = simAnnealOrbSimulator
OBJSsa = simAnnealOrbSimulator.o rnd/sfmt.o rnd/kylesGaussRand.o Toolboxes/orbSimStructures.o Toolboxes/RVtools.o Toolboxes/DItools.o Toolboxes/OrbSimGeneralToolbox.o Toolboxes/SimSettingsObj.o Toolboxes/DataObj.o simAnnealOrbSimFunc.o
	
simAnnealOrbSimulator: $(OBJSsa)
	$(CC) -o simAnnealOrbSimulator $(OBJSsa) $(LDFLAGS) $(LIBS)
#---------------------------------------------------------------
#---------------------------------------------------------------

PROGSmcmc = MCMCorbSimulator
OBJSmcmc = MCMCorbSimulator.o rnd/sfmt.o rnd/kylesGaussRand.o Toolboxes/orbSimStructures.o Toolboxes/RVtools.o Toolboxes/DItools.o Toolboxes/OrbSimGeneralToolbox.o Toolboxes/SimSettingsObj.o Toolboxes/DataObj.o simAnnealOrbSimFunc.o MCMCorbSimFunc.o
	
MCMCorbSimulator: $(OBJSmcmc)
	$(CC) -o MCMCorbSimulator $(OBJSmcmc) $(LDFLAGS) $(LIBS)
#---------------------------------------------------------------	
#---------------------------------------------------------------

PROGSmcmc_l = looped_MCMCorbSimulator
OBJSmcmc_l = looped_MCMCorbSimulator.o rnd/sfmt.o Toolboxes/orbSimStructures.o Toolboxes/RVtools.o Toolboxes/DItools.o Toolboxes/OrbSimGeneralToolbox.o Toolboxes/SimSettingsObj.o Toolboxes/DataObj.o simAnnealOrbSimFunc.o MCMCorbSimFunc.o
	
looped_MCMCorbSimulator: $(OBJSmcmc_l)
	$(CC) -o looped_MCMCorbSimulator $(OBJSmcmc_l) $(LDFLAGS) $(LIBS)
#---------------------------------------------------------------
#
# Update Makefile dependencies
#

-include .Makefile.depend
.PHONY : depend
depend : 
	@echo Rebuilding make dependencies
	$(CC) $(CFLAGS) -MM $(OBJSmc:.o=.cpp) > .Makefile.depend
	$(CC) $(CFLAGS) -MM $(OBJSsa:.o=.cpp) > .Makefile.depend
	$(CC) $(CFLAGS) -MM $(OBJSmcmc:.o=.cpp) > .Makefile.depend
	$(CC) $(CFLAGS) -MM $(OBJSmcmc_l:.o=.cpp) > .Makefile.depend
.PHONY : tags
tags :
	etags *.[ch]pp
#
.PHONY : clean
clean :
	$(RM) $(PROGSmc) gmon.{out,sum} *.o *.s *~ core core.[0-9]*[0-9] TAGS
	$(RM) $(PROGSsa) gmon.{out,sum} *.o *.s *~ core core.[0-9]*[0-9] TAGS
	$(RM) $(PROGSmcmc) gmon.{out,sum} *.o *.s *~ core core.[0-9]*[0-9] TAGS
	$(RM) $(PROGSmcmc_l) gmon.{out,sum} *.o *.s *~ core core.[0-9]*[0-9] TAGS
