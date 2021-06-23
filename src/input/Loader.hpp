#ifndef Loader_hpp
#define Loader_hpp

#include <Python.h>
#include <mpi.h>
#include <stdio.h>
#include <stdexcept>
#include <string>
#include <vector>
#include <math.h>
#include "../misc/Logger.hpp"
#include "../misc/Misc.hpp"


enum RunType{
    SCRATCH,
    RESTART
};

enum FieldBCtype{
    DAMPING,
    PERIODIC
};


enum ParticleBCtype{
    OUTFLOW_BC,
    PERIODIC_BC,
    REFLECT_BC
};


class Loader{
    
private:
    Logger logger;
    double timeStep;
    int numOfSpecies;
    int maxTimestepsNum;
    int timestepsNum2Write;
    int ppc;
    std::string outputDir;
    std::string fileNameTemplate;
    
    PyObject *pInstance;
    
    double callPyFloatFunction( PyObject*, const std::string, const std::string );
    long callPyLongFunction( PyObject*, const std::string, const std::string );
    std::string callPyStringFunction( PyObject*, const std::string, const std::string );
    
    double callPyFloatFunctionWith3args( PyObject*, const std::string,
                                        const std::string, double, double, double );
    
    PyObject * getPyMethod( PyObject*, const std::string, const std::string );
    
    void initMPIcoordinatesOfDomains( int, int[3] );

    
public:
    int runType;
    std::string inputfile;
    int resolution[3];
    int BCtype[3];
    int partclBCtype[3];
    int totPixelsPerBoxSide[3];
    int offsetInPixels[3];
    int mpiDomains[3];
    int mpiCoords[3];
    const int SIMULATION_SIZE = 3;
    int dim;
    double boxCoordinates[3][2];
    
    double dampingBoundaryWidth[3][2];
    
    double boxSizes[3];
    double spatialSteps[3];
    
    double resistivity;

    double minimumDens2ResolvePPC;
    
    int smoothStride;
    double electronmass;
    double relaxFactor;
    
    //stability parameters
    double cellBreakdownEfieldFactor;
    double criticalPressure;
    
    //laser staff
    int numOfSpots = 0;
    int prtclType2Load;
    double prtclTemp2Load;
    double pressureIncreaseRate;
    int laserPulseDuration_tsnum;
    double loadedEnergyPerStep = 0.0;
    
    //MPI staff
    std::vector<int> neighbors2Send;//27
    std::vector<int> neighbors2Recv;//27
    
    
    Loader();
    ~Loader();
    void load();
    double getElectronPressure(double,double,double);
    
    double getTargetIonDensityProfile(double,double,double);
    double getElectronPressureProfile(double,double,double);
    
    double getDensity(double,double,double,int);
    std::vector<double> getVelocity(double,double,double,int);
    std::vector<double> getBfield(double,double,double);
   
    
    double getTimeStep();
    int getMaxTimestepsNum();
    int getNumberOfSpecies();
    int getTimestepsNum2Write();
    
    double getPPC4species(int);
    double getMass4species(int);
    double getCharge4species(int);
    int getIfSpeciesFrozen(int);
    
    std::string getOutputDir() const;
    std::string getFilenameTemplate();
    
    PyObject * getPythonClassInstance(std::string className);
};
#endif /* Loader_hpp */
