
________________________________________________
--------O--------O--A-------A-------55555--2222    
-------O-O------O--A------A-A------5-------------2  
------O---O-----O-A------A---A-----5-------------2  
-----OO-OO----OA------AA--AA----5555----2222    
----O-------O---O-A----A-------A---------5--2----    
---O---------O--O--A--A---------A--55555--22222  
----Arbitrary---Kinetic--Algorithm--------------
________________________________________________

# AKA
 Arbitrary Kinetic Algorithm 
 for collisionless plasma modeling.

 stack: C++11, MPI, HDF5, python2-3

_______________________
#       HOWTO
_______________________
1. before 'make' need to set in makefile    
    HDF5_PATH= path to hdf5 lib    
    MPI_PATH= path to mpi lib     
    PYTHON_INC= path to python include   
    PYTHON_LIB= path to python lib   

2. for running default example from src/input/Initializer.py    
    mpirun -n 2 aka.exe   

3. normally need to set input file path containing Initializer.py     
    mpirun -n 2 aka.exe PATH/TO/PYTHON/INPUT/FILE   

4. before running need to create output folder and set in Initializer.py    

5. for visualization use python notebook in folder NOTEBOOK/    

6. to run it in supercomputer, e.g. Niagara
    `source compile_aka52_niagara.sh`
    edit the `makefile` accordingly, like step 1
    then `make`

___________________________________
#     MPI and HDF5 installation
___________________________________

download openmpi v 4.0.5 from    
download.open-mpi.org/release/open-mpi/v4.0/openmpi-4.0.5.tar.gz   

extract in folder /TEMP/FLD4INSTALLATION/, go to the folder   

run in terminal:    
1. ./configure --prefix=/FLD2LIBS/openmpi CFLAGS=-m64 CXXFLAGS=-m64 FFLAGS=-m64 FCFLAGS=-m64   
// need to specify 64-bit compilation via additional flag that should be passed to ALL compilers   
2. make   
3. make install   

download hdf5 v1.10.5 from    
http://hdfgroup.org/package/hdf5-1-10-5-tar-gz/?wpdmdl=13571    

run in terminal:   
1. ./configure --prefix=/FLD2LIBS/hdf5 --enable-parallel CC=/FLD2LIBS/openmpi/bin/mpicc   
2. make   
3. make install   

_______________________
# physics included:
_______________________   
* space and time quantities are normalized on:  
  - ion inertia length d0 = c/omega_pi  
  - inverse ion gyrofrequency 1/Omega_ci  

* electro-magnetic fields are calculated on two staggered grid (G1 and G2)  
  using predictor-corrector scheme (see EleMagManager.cpp)  

  E = -[VixB] + [JxB]/n - divPe/n - eta J   

  B' = -rotE  

  J = rotB  

* ions are described in a kinetic way,   
  dynamics is solved using first order interpolation of em fields  

* electrons are described in a fluid way by:   
  - density (equals to ion density ne=ni=n)   
  - bulk velocity Ve = Vi - J/n  
  - six-component pressure tensor Pij  

* six-component pressure tensor P is integrated in time   
  using subcycling explicit scheme, 
  for implicit one need to build with flag -DIMPLICIT_PRESSURE

  P' = - Ve.∇P - P∇.Ve - P.∇Ve - (P.∇Ve)^T - q/m [ PxB + (PxB)^T ]   

  where Ve - electron flow velocity  
        q - electron charge    
        m - electron mass   
        B - magnetic field   
 
* laser effects are imitated by ablation operator including  
  heat operator and particle creation operator   
  (see LaserMockManager.cpp)  

* ablation operator works in localized area, called focal spot  

* heat operator provides electron pressure increasing in the focal spot  

* particle creation operator sustains constant target density   



_______________________
#     FEATURES
_______________________
1. AKA is 3D, parallel (MPI), multispecies code 

2. BC type for em fields and hydro quantities: 
   1 - periodic (see GridManager.cpp)
   0 - damping layer (see EleMagManager.cpp and ClosureManager.cpp)

3. BC type for particle properties (see BoundaryManager.cpp):  
   1 - periodic 
   0 - outflow // reaching border particle leaves domain forever

4. for small scale dissipation use resistivity (eta) parameter

5. for pressure tensor integration need to set:   
   * electron mass   
   * relaxation factor for izotropization operator (see ClosureManager.cpp)     
   * smooth stride for pressure tensor smoothing    

6. use 'make FLAGS=-DLOG' to  see some logs    
   to set debug log level see Logger.hpp   
   for extra logging use -DHEAVYLOG    

7. for each particle type have to specify:   
   * mass   
   * charge   
   * if particles are frozen in space (skip pusher phase)  
   * particles per cell number   

8. if number of laser focal spots is more than zero,      
   auxiliary particle type is reserved for the loading fraction 

_______________________
#   TROUBLESHOUTING:
_______________________
1.  for OS X: ....python2.7/pyport.h:731:29:   
    note: expanded from macro 'toupper'    
    solution:     
    add #ifndef __cplusplus... line 699 + #endif    
    in /pyport.h line 722    


