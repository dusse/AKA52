HDF5_PATH=$(MODULE_HDF5_PREFIX)
MPI_PATH=$(MODULE_OPENMPI_PREFIX)
PYTHON27_INC=$(SCINET_PYTHON_ROOT)/include/python2.7
PYTHON27_LIB=$(SCINET_PYTHON_ROOT)/lib


LIBS=-lpython2.7 -lhdf5
INCLUDES=-I$(PYTHON27_INC) -I$(HDF5_PATH)/include/ -I$(MPI_PATH)/include/

DSRC = ./src
DEXE = ./

LD_LIBRARY_PATH=$(HDF5_PATH)/lib/:$(PYTHON27_LIB)


export LIBRARY_PATH=$LIBRARY_PATH:$(LD_LIBRARY_PATH)

CXX = $(MPI_PATH)/bin/mpicxx
CXXFLAGS  = -Wall -c -std=c++11 -Wno-sign-compare -Wno-unused-variable

_SRCS =  $(DSRC)/core/SimulationManager.cpp \
               $(DSRC)/grid/GridManager.cpp \
               $(DSRC)/grid/boundary/BoundaryManager.cpp \
               $(DSRC)/input/Loader.cpp \
	       $(DSRC)/particles/Particle.cpp \
               $(DSRC)/misc/Logger.cpp \
               $(DSRC)/misc/Misc.cpp \
               $(DSRC)/output/Writer.cpp \
               $(DSRC)/physics/pusher/Pusher.cpp \
               $(DSRC)/physics/hydro/HydroManager.cpp \
               $(DSRC)/physics/electro-magnetic/EleMagManager.cpp \
               $(DSRC)/physics/pressure-closure/ClosureManager.cpp \
               $(DSRC)/physics/laser/LaserMockManager.cpp \
               $(DSRC)/physics/collisions/IonIonCollisionManager.cpp \
               $(DSRC)/common/variables/VectorVar.cpp \
               $(DSRC)/solvers/Solver.cpp \
               $(DSRC)/solvers/ModelInitializer.cpp \
               $(DSRC)/AKA.cpp \

_OBJS            = $(_SRCS:.cpp=.o)

_EXEN            = $(DEXE)/aka.exe

all : $(_EXEN)


$(_EXEN) : $(_OBJS)
	@echo 'Building target: $@'
	$(CXX) -o $@ $^  $(LIBS)
	@echo 'Finished building target: $@'
	@echo ' '

%.o : %.cpp
	$(CXX) $(INCLUDES) -o $@ -c $< $(CXXFLAGS) $(FLAGS)


clean :
	rm -f $(_OBJS)


