 #***********************************************************************
 # $Id: Makefile,v 1.1.1.1 2003/12/11 08:18:16 tprosser Exp $
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
#CPP = icc 
#CPP  = /home/tprosser/gcc/bin/g++
CPP = g++  
CFLAGS = -fabi-version=0 -Wno-long-double  
RM  = rm
RMFLAGS = -fr
LIB =-framework vecLib  -lfftw3 -lm 
#LIBS = -L/home/tprosser/libraries/LAPACK -llapack -lblas -lm -lg2c
#LIBS2 = $(LIB) -L/home/tprosser/gcc/lib -lg2c -L/home/tprosser/libraries/fftw/lib -lfftw3
#Blitz should be installed with matching compiler option corresponding to the specific compiler.
#LIBSS = $(LIBS) -L/home/tprosser/libraries/blitz/lib/ -lblitz
#LIBSS = $(LIBS) -L$(HOME)/blitz-0.5beta3/lib -lblitz
#INCLDIR = /home/tprosser/gcc/include/
BOOSTDIR =/usr/local/include/boost/      
ITLDIR = /usr/local/include/itl      
MTLDIR = /usr/local/include/mtl            
IETLDIR =/usr/local/include/ietl       
FFTWDIR=/usr/local/include/      
#BLITZDIR = /home/tprosser/libraries/blitz/include/
#BLITZDIR = $(HOME)/blitz-0.5beta3

.SUFFIXES:      .cpp      .o

.cpp.o:
	@$(ECHO) Making $@ from $<
	@$(CPP) -c $(CFLAGS) -I$(BOOSTDIR) -I$(ITLDIR) -I$(IETLDIR) -I$(BLITZDIR) -I$(MTLDIR) -I$(FFTWDIR) $<

all:   BEC_Groundstate

BEC_Groundstate:BEC_Groundstate.o
	@$(ECHO)        Linking
	@$(CPP) -o $@ BEC_Groundstate.o $(LIB)


BEC_try:BEC_try.o                
	@$(ECHO)        Linking
	@$(CPP) -o $@ BEC_try.o $(LIB)


BEC_2try:BEC_2try.o                
	@$(ECHO)        Linking
	@$(CPP) -o $@ BEC_2try.o $(LIB)


BEC_3try:BEC_3try.o                
	@$(ECHO)        Linking
	@$(CPP) -o $@ BEC_3try.o $(LIB)


BEC_harm_lattice_nokick:BEC_harm_lattice_nokick.o
	@$(ECHO)        Linking
	@$(CPP) -o $@ BEC_harm_lattice_nokick.o $(LIB)


BEC_harm_lattice_nokick_new:BEC_harm_lattice_nokick_new.o
	@$(ECHO)        Linking
	@$(CPP) -o $@ BEC_harm_lattice_nokick_new.o $(LIB)

BEC_harm_lattice_kick:BEC_harm_lattice_kick.o
	@$(ECHO)        Linking
	@$(CPP) -o $@ BEC_harm_lattice_kick.o $(LIB)

BEC_harm_lattice_kick1:BEC_harm_lattice_kick1.o
	@$(ECHO)        Linking
	@$(CPP) -o $@ BEC_harm_lattice_kick1.o $(LIB)

BEC_kick_mom:BEC_kick_mom.o
	@$(ECHO)        Linking
	@$(CPP) -o $@ BEC_kick_mom.o $(LIB)

BEC_Energy:BEC_Energy.o
	@$(ECHO)        Linking
	@$(CPP) -o $@ BEC_Energy.o $(LIB)


Mean_BEC_momentum:Mean_BEC_momentum.o
	@$(ECHO)        Linking
	@$(CPP) -o $@ Mean_BEC_momentum.o $(LIB)


Mean_BEC_momentum_dist:Mean_BEC_momentum_dist.o
	@$(ECHO)        Linking
	@$(CPP) -o $@ Mean_BEC_momentum_dist.o $(LIB)

Mean_BEC_momentum_dist_notrap:Mean_BEC_momentum_dist_notrap.o
	@$(ECHO)        Linking
	@$(CPP) -o $@ Mean_BEC_momentum_dist_notrap.o $(LIB)

Mean_BEC_momentum_dist_nothing:Mean_BEC_momentum_dist_nothing.o
	@$(ECHO)        Linking
	@$(CPP) -o $@ Mean_BEC_momentum_dist_nothing.o $(LIB)

legendre:legendre.o                            
	@$(ECHO)        Linking
	@$(CPP) -o $@ legendre.o $(LIB)

BEC_Groundstate_imt:BEC_Groundstate_imt.o
	@$(ECHO)        Linking
	@$(CPP) -o $@ BEC_Groundstate_imt.o $(LIB)

clean:
	@$(RM) $(RMFLAGS) *.o ti_files ring
