#########################################################
###        Part A: change with your environment       ###
#########################################################
###1. system platform
PLATFORM = soilblast
#PLATFORM = spirit
#PLATFORM = garnet
#PLATFORM = janus
#PLATFORM = lightning
#PLATFORM = excalibur

###2. compiler combinations (c++ and mpi)
COMPILERS = gcc-openmpi
#COMPILERS = intel-intel
#COMPILERS = intel-sgi
#COMPILERS = intel-cray

###3. boost location
ifeq ($(PLATFORM), soilblast)
 BOOST_ROOT = /usr/local/boost-1.53.0-openmpi-1.6.4-gcc-4.6.4
 #BOOST_ROOT = default
endif

ifeq ($(PLATFORM), spirit)
 ifeq ($(COMPILERS), gcc-openmpi)
  BOOST_ROOT = /work4/projects/openmpi/boost-1.55.0-openmpi-1.7.5-gcc-compilers-4.7.3
 endif
 ifeq ($(COMPILERS), intel-intel)
  #BOOST_ROOT = /work4/projects/openmpi/boost-1.53.0-intelmpi-4.0.3.008
  BOOST_ROOT = /work4/projects/openmpi/boost-1.55.0-intelmpi-4.1.3.048-intel-compilers-13.0.1
 endif
 ifeq ($(COMPILERS), intel-sgi)
  BOOST_ROOT = /work4/projects/openmpi/boost-1.57.0_mpt-2.11_intel-compilers-13.0.1
  #BOOST_ROOT = /home/yanb/local/boost-1.57.0_mpt-2.11_intel-compilers-13.0.1
 endif
endif

ifeq ($(PLATFORM), garnet)
 ifeq ($(COMPILERS), intel-cray)
  BOOST_ROOT = /u/yanb/local/boost-1.55.0_cray-mpich6.0.0_PrgEnv-intel4.2.34
 endif
endif

ifeq ($(PLATFORM), lightning)
 ifeq ($(COMPILERS), intel-cray)
  BOOST_ROOT = /home/yanb/local/boost-1.55.0_cray-mpich6.3.1_PrgEnv-intel5.2.14
 endif
endif

ifeq ($(PLATFORM), excalibur)
 ifeq ($(COMPILERS), intel-cray)
  BOOST_ROOT = /usr/people/yanb/local/boost-1.57.0_cray-mpich7.1.0_PrgEnv-intel5.2.40
 endif
endif

ifeq ($(PLATFORM), janus)
 ifeq ($(COMPILERS), gcc-openmpi)
  #BOOST_ROOT = /home/yanb/local/boost-1.55.0_openmpi-1.7.3_gcc-4.8.2_ib
  BOOST_ROOT = /home/yanb/local/boost-1.55.0_openmpi-1.8.2_gcc-4.9.1
 endif
 ifeq ($(COMPILERS), intel-intel)
  BOOST_ROOT = /home/yanb/local/boost-1.55.0_impi-13.0.0_intel-13.0.0
 endif
endif

#########################################################
###         Part B: usually do not change             ###
#########################################################
###1. MPI wrapper, OpenMP and optimize/debug
# MPI wrapper
ifeq ($(COMPILERS), gcc-openmpi)
 MPICXX = mpic++
endif

ifeq ($(COMPILERS), intel-intel)
 MPICXX = mpiicpc
endif

ifeq ($(COMPILERS), intel-cray)
 MPICXX = CC
endif

ifeq ($(COMPILERS), intel-sgi)
 MPICXX = mpicxx
endif

# OpenMP
ifeq ($(COMPILERS), gcc-openmpi)
 OPENMP = -fopenmp
endif

ifeq ($(COMPILERS), intel-intel)
 OPENMP = -openmp -mt_mpi
endif

ifeq ($(COMPILERS), intel-cray)
 OPENMP = -openmp
endif

ifeq ($(COMPILERS), intel-sgi)
 OPENMP = -openmp -mt_mpi
endif

# optimize or debug
OPTIMIZE = -O3 -DNDEBUG 
#OPTIMIZE = -Wall -g
#gcc: full level of IEEE 754 compliance
#OPTIMIZE = -O3 -DNDEBUG -frounding-math -fsignaling-nans
#intel: -fp-model source also implies keyword precise; and it is nearly as twice slow as default -fp-model fast=1
#OPTIMIZE = -O3 -DNDEBUG -fp-model source

# GCC quadmath, default off
#QUADMATH = -DQUADMATH -L/usr/local/gcc-4.6.2/lib64 -lquadmath

# CXXFLAGS
CXXFLAGS = $(OPTIMIZE) $(OPENMP) $(QUADMATH)

###2. BOOST details
ifeq ($(PLATFORM), soilblast)
  ifeq ($(BOOST_ROOT), default)
   BOOST_LIB_DIR=-L/usr/lib64/openmpi/lib:/usr/lib64/compat-openmpi/lib:/usr/lib64
   BOOST_LIBS=-lboost_mpi -lboost_serialization -lboost_timer -lboost_system -lboost_chrono
   BOOST_RUN_LIB_DIR=-Wl,-rpath=/usr/lib64/openmpi/lib:/usr/lib64/compat-openmpi/lib:/usr/lib64
  else
   BOOST_INCLUDE=-I$(BOOST_ROOT)/include
   BOOST_LIB_DIR=-L$(BOOST_ROOT)/lib
   BOOST_LIBS=-lboost_mpi -lboost_serialization -lboost_timer -lboost_system -lboost_chrono
   BOOST_RUN_LIB_DIR=-Wl,-rpath=$(BOOST_ROOT)/lib
  endif
else
# for spirit, janus, garnet
 BOOST_INCLUDE=-I$(BOOST_ROOT)/include
 BOOST_LIB_DIR=-L$(BOOST_ROOT)/lib
 BOOST_LIBS=-lboost_mpi -lboost_serialization -lboost_timer -lboost_system -lboost_chrono
 BOOST_RUN_LIB_DIR=-Wl,-rpath=$(BOOST_ROOT)/lib
endif

INCLUDE=$(BOOST_INCLUDE)
LIB_DIR=$(BOOST_LIB_DIR)
LIBS=$(BOOST_LIBS)
RUN_LIB_DIR=$(BOOST_RUN_LIB_DIR)

###3. makefile
SOURCES = $(wildcard *.cpp)
OBJECTS = $(SOURCES:.cpp=.o)
#OBJECTS = $(patsubst %.cpp, %.o, $(SOURCES))
EXECUTABLE = paraEllip3d

.PHONY: all tar clean

all: $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(MPICXX) -o $@ $(CXXFLAGS) $(OBJECTS) $(INCLUDE) $(LIB_DIR) $(LIBS) $(RUN_LIB_DIR)

%.o: %.cpp
	$(MPICXX) -c $< -o $@ $(CXXFLAGS) $(INCLUDE) $(LIB_DIR) $(LIBS)

sinclude $(SOURCES:.cpp=.d)

%.d: %.cpp
	$(MPICXX) -MM $< > $@.$$$$; \
	sed 's,\($*\)\.o[ :]*,\1.o $@ : ,g' < $@.$$$$ > $@; \
	rm -f $@.$$$$

tar:
	tar -cvf $(EXECUTABLE).tar *.h *.cpp makefile* readme

clean:
	-rm -f *.o *.d  *~ *.tar $(EXECUTABLE)

# sinclude is always resolved even if make tar/clean

