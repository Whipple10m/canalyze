#*****************************************************************************
# Makefile for Purdue Whipple Data Analysis Software
# Written by R.W. Lessard
#*****************************************************************************
# Make changes here only
#
SHELL = /bin/csh
#
# Change your compiler here (e.g. egcs, mpicc)
# If you don't want MPI remove -DMPI and set CC to your compiler
CC = gcc
CFLAGS = -O 
FC = g77
FFLAGS = -O
# MPI run command 
# If no MPI
MPI_RUN = 
#MPI_RUN = "mpirun -np 1"
#
LDFLAGS =
#
# Define your machine name
# options SUN,DEC,LINUX,VMS
ARCH = LINUX
#
# Define place where data is stored (usage like PATH variable, : delimited)
#
DATAPATH=.:/cdrom:/data
#
# Define place where database is stored
#
DBPATH=.
#
# Path to install binaries
#
INSTALL_DIR=/home/finley/bin4.0
#
# Define path and names of pulsar and binary catalogs
#
PSRCAT = "/home/finley/bin/psrjcat.dat"
CRABCAT = "/home/finley/bin/crabephem.dat"
BINCAT = "/home/finley/bin/binarycat.dat"
#*****************************************************************************
#
# Changes below only required to define path to system and external libraries 
#
# Change the path to your libraries
# X-windows
LIBX = -L/usr/X11R6/lib
LIB_PREFIX = $(HOME)/src
# numerical recipes
LIBREC = -L$(LIB_PREFIX)/recipes 
INCREC = -I$(LIB_PREFIX)/recipes
# granite data format (required only for fz2hdf)
LIBGDF = -L$(LIB_PREFIX)/gdf/vdev
GDF = -lgdf 
# cernlib (required only for fz2hdf)
LIBCERN = -L$(LIB_PREFIX)/cern/99/lib
CERN = -lmathlib -lpacklib -lkernlib
# starlink - SLA
LIBSTAR = -L$(LIB_PREFIX)/star-i686/lib
# HDF
LIBHDF = -L$(LIB_PREFIX)/4.1r3_linux/lib
INCHDF = -I$(LIB_PREFIX)/4.1r3_linux/include
HDF = -ldf -ljpeg -lz
# MATLAB
LIBMAT = -L/usr/local/matlabr12/extern/lib/glnx86
INCMAT = -I/usr/local/matlabr12/extern/include
MAT = -lmat -lmx -lmi -lut
# Barycentering code by J.P. Finley
LIBBARY  = -L$(LIB_PREFIX)/bary
BARY = -lbary
#
# System dependent libraries
LIBS = -L. -lfio -lf77math -lf90math -lg2c
#*****************************************************************************
# your done - changes below unlikely 

cpeds: cpeds.o cpeds.c whipple.o
	$(CC) cpeds.o whipple.o $(CFLAGS) $(LDFLAGS) -o cpeds \
	$(LIBHDF) $(HDF) \
	$(LIBSTAR) -lsla \
	$(LIBS) -lm

cpeds.o: cpeds.c
	$(CC) cpeds.c $(CFLAGS) -c -D$(ARCH) \
	$(INCHDF)

cn2gains: cn2gains.o cn2gains.c whipple.o
	$(CC) cn2gains.o whipple.o $(CFLAGS) $(LDFLAGS) -o cn2gains \
	$(LIBHDF) $(HDF) \
	$(LIBSTAR) -lsla \
	$(LIBS) -lm

cn2gains.o: cn2gains.c
	$(CC) cn2gains.c $(CFLAGS) -c -D$(ARCH) \
	$(INCHDF)

cparam: cparam.o cparam.c whipple.o sub_macro_wv_analysis.o
	$(CC) cparam.o whipple.o sub_macro_wv_analysis.o \
	$(CFLAGS) $(LDFLAGS) -o cparam \
	$(LIBHDF) $(HDF) \
	$(LIBSTAR) -lsla \
	$(LIBREC) -lrec \
	$(LIBS) -lm

cparam.o: cparam.c
	$(CC) cparam.c $(CFLAGS) -c -D$(ARCH) \
	$(INCHDF) $(INCREC)

sub_macro_wv_analysis.o: sub_macro_wv_analysis.f
	$(FC) sub_macro_wv_analysis.f $(FFLAGS) -c \

mcut: mcut.o mcut.c mbook.o whipple.o 
	$(CC) mcut.o mbook.o whipple.o $(CFLAGS) $(LDFLAGS) -o mcut \
	$(LIBHDF) $(HDF) \
	$(LIBMAT) $(MAT) \
	$(LIBSTAR) -lsla \
	$(LIBBARY) $(BARY) \
	$(LIBS) -lm

mcut.o: mcut.c
	$(CC) mcut.c -c $(CFLAGS) -D$(ARCH) \
	-DHAVE_MATLAB $(INCMAT) \
	$(INCHDF)

mbook.o: mbook.c
	$(CC) mbook.c -c $(CFLAGS) \
	$(INCMAT)

find_ephemeris: find_ephemeris.o find_ephemeris.c
	$(CC) find_ephemeris.o $(CFLAGS) $(LDFLAGS) -o find_ephemeris \
	$(LIBSTAR) -lsla $(LIBS)

find_ephemeris.o: find_ephemeris.c
	$(CC) find_ephemeris.c $(CFLAGS) -c -D$(ARCH) \
	-DPULSAR_CATALOG='$(PSRCAT)' -DCRAB_CATALOG='$(CRABCAT)'

find_orbit: find_orbit.o find_orbit.c
	$(CC) find_orbit.o $(CFLAGS) $(LDFLAGS) -o find_orbit \
	$(LIBSTAR) -lsla $(LIBS)

find_orbit.o: find_orbit.c
	$(CC) find_orbit.c $(CFLAGS) -c -D$(ARCH) \
	-DBINARY_CATALOG='$(BINCAT)'

fz2hdf: fz2hdf.o mbook.o fz2hdf.c
	$(CC) fz2hdf.o mbook.o $(CFLAGS) $(LDFLAGS) -o fz2hdf \
	$(LIBGDF) $(LIBCERN) $(LIBHDF) $(LIBMAT) \
	$(LIBSTAR) -lsla $(MAT) $(GDF) $(CERN) $(HDF) \
	$(LIBS) -lm

fz2hdf.o: fz2hdf.c
	$(CC) fz2hdf.c $(CFLAGS) -D$(ARCH) -DHAVE_MATLAB -c \
	$(INCMAT) $(INCHDF) $(INCMBOOK)

red2hdf: red2hdf.o fileio.o vtoi.o whipple.o red2hdf.c
	$(CC) red2hdf.o fileio.o vtoi.o whipple.o $(CFLAGS) \
	$(LDFLAGS) -o red2hdf \
	$(LIBHDF) $(HDF) \
	$(LIBS) -lm

red2hdf.o: red2hdf.c
	$(CC) red2hdf.c $(CFLAGS) -D$(ARCH) -c $(INCHDF) 

fileio.o: fileio.c
	$(CC) fileio.c $(CFLAGS) -D$(ARCH) -c

vtoi.o: vtoi.c
	$(CC) vtoi.c $(CFLAGS) -D$(ARCH) -c

canalyze.config:
	echo DEFAULTDATAPATH=$(DATAPATH) > canalyze.config ; \
	echo DBPATH=$(DBPATH) >> canalyze.config ; \
	echo BINPATH=$(INSTALL_DIR) >> canalyze.config
	echo MPICOMM='$(MPI_RUN)' >> canalyze.config

all: red2hdf fz2hdf cpeds cn2gains cparam mcut find_ephemeris find_orbit canalyze.config

clean: 
	rm -f *.o core ; \
	rm -f red2hdf fz2hdf cpeds cn2gains cparam mcut; \
	rm -f find_ephemeris find_orbit canalyze.config

install: 
	cp red2hdf $(INSTALL_DIR); \
	cp fz2hdf $(INSTALL_DIR); \
	cp cpeds $(INSTALL_DIR); \
	cp cn2gains $(INSTALL_DIR); \
	cp cparam $(INSTALL_DIR); \
	cp mcut $(INSTALL_DIR); \
	cp find_ephemeris $(INSTALL_DIR); \
	cp find_orbit $(INSTALL_DIR); \
	cp canalyze $(INSTALL_DIR); \
	cp canalyze.config $(INSTALL_DIR); \
	cp *.m $(INSTALL_DIR); \
	cp *.mat $(INSTALL_DIR); \
	cp *.dat $(INSTALL_DIR)
