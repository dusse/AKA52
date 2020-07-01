
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

 stack: C++11, MPI, HDF5, python2.7

_______________________
#       HOWTO
_______________________
1. before 'make' need to set in makefile  
    HDF5_PATH= path to hdf5 lib (last well used 1.10.5)  
    MPI_PATH= path to mpi lib (last well used openmpi 9.0.0)  
    PYTHON27_INC= path to python2.7 include  
    PYTHON27_LIB= path to python2.7 lib  

2. for running default example (harris 2D, GEM challenge) from src/input/Initializer.py  
    mpirun -n 2 aka.exe

3. normally need to set input file path containing Initializer.py  
    mpirun -n 2 aka.exe PATH/TO/PYTHON/INPUT/FILE

4. before running need to create output folder and set in Initializer.py  

5. for visualization use python notebook in folder NOTEBOOK/  

_______________________
# physics included:
_______________________
* space and time quantities are normalized on:  
  - ion inertia length d0 = c/omega_pi  
  - inverse ion gyrofrequency 1/Omega_ci  

* electro-magnetic fields are calculated on two staggered grid (G1 and G2)  
  using predictor-corrector scheme (see EleMagManager.cpp)  

  E = -[VixB] + [JxB]/n - divPe/n - eta ΔJ   

  B' = -rotE  

  J = rotB  

* ions are described in a kinetic way,   
  dynamics is solved using first order interpolation of em fields  

* electrons are described in a fluid way by:   
  - density (equals to ion density ne=ni=n)   
  - bulk velocity Ve = Vi - J/n  
  - six-component pressure tensor Pij  

* six-component pressure tensor P is integrated in time   
  using subcycling explicit scheme  

  P' = - Ve.∇P - P∇.Ve - P.∇Ve - (P.∇Ve)^T - q/m [ PxB + (PxB)^T ]   

  where Ve - electron flow velocity  
        q - electron charge  
        m - electron mass  
        B - magnetic field  
 
* laser effects are imitated by ablation operator including  
  heat operator and particle creation operator   
  (see LaserMockManager.cpp)  

* ablation operator works in localized area, called focal spot  

* heat operator provides ions heating and pressure increasing in the focal spot  

* particle creation operator sustains constant target density   



_______________________
#     FEATURES
_______________________
1. AKA is 3D, parallel (MPI), multispecies code 

2. BC type for em fields and hydro quantities: 1 - periodic, 0 - ideal (see GridManager.cpp)

3. BC type for particle properties: 2 - reflect 1 - periodic 0 - outflow (see BoundaryManager.cpp)
   *outflow BC - reaching border particle leaves domain forever

4. for small scale dissipation use hyperviscosity (eta) parameter

5. for pressure tensor integration need to set:
   * electron mass
   * relaxation factor for izotropization operator (see ClosureManager.cpp)
   * smooth stride for pressure tensor smoothing

6. use 'make FLAGS=-DLOG' to set debug log level in Logger.hpp


_______________________
#        TODO:
_______________________
1. [2] fix restart file writing for more than 2 cores
2. [0] include collisions



_______________________
#   TROUBLESHOUTING:
_______________________
1.  for OS X: ....python2.7/pyport.h:731:29:  note: expanded from macro 'toupper'
    solution: add #ifndef __cplusplus...line 699 + #endif in /pyport.h line 722







