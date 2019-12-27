.POSIX:
#----------------------------------------
#
#  file created 98jun18, cca, (cleve.ashcraft@boeing.com)
#  based on work by clay breshears (clay@turing.wes.hpc.mil)
#  (much appreciated)
#
#----------------------------------------
#
#  place your favorite compiler here
#
# for solaris
#
 CC = gcc
#  CC = /usr/lang-4.0/bin/cc
#
# for sgi
#
# CC = cc
#
# for hp
#
# CC = /opt/mpi/bin/mpicc
#
#----------------------------------------
#
#  set the compiler flags
# 
# OPTLEVEL = 
# OPTLEVEL = -g -v
  OPTLEVEL = -O
# OPTLEVEL = -xO5 -v
# OPTLEVEL = -O3
# OPTLEVEL = -O4
# CFLAGS = -Wall -g
# CFLAGS = -Wall -pg
# CFLAGS = $(OPTLEVEL) -D_POSIX_C_SOURCE=199506L
  CFLAGS = $(OPTLEVEL)
# CFLAGS = -Wall $(OPTLEVEL)
#
#----------------------------------------
#
#  set any load flags
#
#  LDFLAGS = -Wl,+parallel -Wl,+tm,spp2000 # for hp exemplar
   LDFLAGS =
#
#---------------------------------------------------------------------
#
#  set any thread libraries
#
# THREAD_LIBS =
# THREAD_LIBS = -D_REENTRANT=199506L -lpthread
  THREAD_LIBS = -D_POSIX_C_SOURCE=199506L -lpthread
# THREAD_LIBS = -lpthread
#
#---------------------------------------------------------------------
#
#  set the purify environment (a memory monitoring tool)
#
  PURIFY = 
# PURIFY = /usr/local/purify-4.0.1/purify
#
#  purify wouldn't work with the newest version of the gcc library,
#  so we had to force loading the old version
#
  PURIFY_GCC_VERSION = 
# PURIFY_GCC_VERSION = -V 2.7.2
#
#---------------------------------------------------------------------
#
#  set the archive flags
#
  AR = ar
  ARFLAGS = rv
#
#---------------------------------------------------------------------
#
#  set the ranlib environment
#  (if ranlib is not needed, we echo the library name)
#
# RANLIB = ranlib
  RANLIB = echo
#
#---------------------------------------------------------------------
#
#  set suffix rule *.c --> *.o
#
.c.o : 
	$(PURIFY) $(CC) -c $(CFLAGS) $<
#
#---------------------------------------------------------------------
#
#  set suffix rule *.c --> *.a
#
.c.a : 
	$(PURIFY) $(CC) -c $(CFLAGS) $<
	$(AR) $(ARFLAGS) $@ $*.o
	rm -f $*.o
#
#---------------------------------------------------------------------
#
#  MPI install library
#
# MPI_INSTALL_DIR = 
  MPI_INSTALL_DIR = /usr/local/mpich-1.0.13
#
#---------------------------------------------------------------------
#
#  MPI library path
#
# for sgi
#
# MPI_LIB_PATH = 
#
# for solaris
#
  MPI_LIB_PATH = -L$(MPI_INSTALL_DIR)/lib/solaris/ch_p4
#
# for hp
#
# MPI_LIB_PATH = 
#
#---------------------------------------------------------------------
#
#  MPI libraries 
#
# for solaris
#
  MPI_LIBS = $(MPI_LIB_PATH) -D_REENTRANT -lmpi -lsocket -lnsl -lthread
#
# for sgi
#
# MPI_LIBS = $(MPI_LIB_PATH) -lmpi -lpthread 
#
# for hp
# MPI_LIBS = -lpthread
# MPI_LIBS = $(MPI_LIB_PATH) -lpthread 
#
#---------------------------------------------------------------------
#
#  MPI include path
#
# MPI_INCLUDE_DIR = 
  MPI_INCLUDE_DIR = -I$(MPI_INSTALL_DIR)/include
#
#---------------------------------------------------------------------
