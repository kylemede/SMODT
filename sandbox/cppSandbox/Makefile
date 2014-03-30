.c.o:
	$(CC) -c $(CFLAGS) $*.c

CC = g++ -std=gnu++11
CFLAGS = -Wall -pedantic -fpermissive
LDFLAGS = -pg
LIBS = -lm -lrt

PROGS = helloWorld

OBJS = rnd/stoc1.o helloWorld.o headerTest.o 

helloWorld: $(OBJS)
	$(CC) -o helloWorld $(OBJS) $(LDFLAGS) $(LIBS)
#
# Update Makefile dependencies
#

PROGSmc = mcONLYorbSimTH_I
OBJSmc = mcONLYorbSimTH_I.o orbSimStructures.o VRtoolbox.o DItoolbox.o OrbSimGeneralToolbox.o

mcONLYorbSimTH_I: $(OBJSmc)
	$(CC) -o mcONLYorbSimTH_I $(OBJSmc) $(LDFLAGS) $(LIBS)
#
# Update Makefile dependencies
#


-include .Makefile.depend
.PHONY : depend
depend : 
	@echo Rebuilding make dependencies
	$(CC) $(CFLAGS) -MM $(OBJS:.o=.cpp) > .Makefile.depend
.PHONY : tags
tags :
	etags *.[ch]pp
#
.PHONY : clean
clean :
	$(RM) $(PROGS) gmon.{out,sum} *.o *.s *~ core core.[0-9]*[0-9] TAGS
