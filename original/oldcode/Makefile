 #***********************************************************************
 # $Id: Makefile,v 1.1.1.1 2003/12/11 08:18:17 tprosser Exp $
 #
 # Makefile, for building the example programs.
 #
 # Copyright (C) 2001-2002 by Prakash Dayal <prakash@comp-phys.org>
 #                            Matthias Troyer <troyer@comp-phys.org>
 #
 # This program is free software; you can redistribute it and/or
 # modify it under the terms of the GNU General Public License
 # as published by the Free Software Foundation; either version 2
 # of the License, or (at your option) any later version.
 #
 # This program is distributed in the hope that it will be useful,
 # but WITHOUT ANY WARRANTY; without even the implied warranty of
 # MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 # GNU General Public License for more details.
 #
 #*************************************************************************/

ECHO = echo
#CPP = g++-3.3 -fabi-version=0
CPP  = icc
#CPP = g++  
CFLAGS =-DNDEBUG 
RM  = rm
RMFLAGS = -fr
LIB = -lm -llapack -lblas
LIBS= $(LIB) -lf2c
#Blitz should be installed with matching compiler option corresponding to the specific compiler.
#LIBSS = $(LIBS) -L/home/tprosser/libraries/blitz/lib/ -lblitz
#LIBSS = $(LIBS) -L$(HOME)/blitz-0.5beta3/lib -lblitz
#INCLDIR = /home/tprosser/gcc/include/
BOOSTDIR = /home/tprosser/libraries/boost/
ITLDIR = /home/tprosser/libraries/itl/
MTLDIR = /home/tprosser/libraries/include/
IETLDIR =/home/tprosser/ietl/include/
#BLITZDIR = /home/tprosser/libraries/blitz/include/
#BLITZDIR = $(HOME)/blitz-0.5beta3

.SUFFIXES:      .cpp      .o

.cpp.o:
	@$(ECHO) Making $@ from $<
	$(CPP) -c $(CFLAGS) -I$(BOOSTDIR) -I$(ITLDIR) -I$(IETLDIR) -I$(BLITZDIR) -I$(MTLDIR) $<

all:   BEC_Groundstate

BEC_Groundstate: BEC_Groundstate.o
	@$(ECHO)        Linking
	@$(CPP) -o $@ BEC_Groundstate.o $(LIBS)



clean:
	@$(RM) $(RMFLAGS) *.o ti_files ring
