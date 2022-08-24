#ifndef Loader_hpp
#define Loader_hpp

#include <Python.h>
#include <mpi.h>
#include <stdio.h>
#include <fstream>
#include <stdexcept>
#include <string>
#include <vector>
#include <set>
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

#define METHOD_OK   0
#define METHOD_FAIL 1

class Loader{
    
private:
    Logger logger;
    double timeStep;
    int numOfSpecies;
    int maxTimestepsNum;
    int timestepsNum2Write;
    std::string outputDir;
    std::string fileNameTemplate;
    std::set<std::string> checkedMethods;
    std::set<std::string> failMethods;
    
    double  collisionFrequencyFactor;
    double defaultCoulombLogarithm;

    PyObject *pInstance;
    
    double callPyFloatFunction( PyObject*, const std::string, const std::string );
    long callPyLongFunction( PyObject*, const std::string, const std::string );
    std::string callPyStringFunction( PyObject*, const std::string, const std::string );
    
    double callPyFloatFunctionWith3args( PyObject*, const std::string,
                                        const std::string, double, double, double );
    
    PyObject * getPyMethod( PyObject*, const std::string, const std::string );
    
    void initMPIcoordinatesOfDomains( int, int[3] );
    int checkMethodExistence(const std::string, const std::string);


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

    int useIsothermalClosure = 0;
    double electronTemperature = 0.0;
    
    int smoothStride;
    double electronmass;
    double relaxFactor;
    
    //stability parameters
    double cellBreakdownEfieldFactor;
    double criticalPressure;
    
    //laser staff
    int numOfSpots = 0;
    int prtclType2Load;
    int laserPulseDuration_tsnum;
    double loadedEnergyPerStep = 0.0;
    
    //MPI staff
    std::vector<int> neighbors2Send;//27
    std::vector<int> neighbors2Recv;//27
    
    
    Loader();
    ~Loader();
    void load();
    double getElectronPressureXX(double,double,double);
    double getElectronPressureYY(double,double,double);
    double getElectronPressureZZ(double,double,double);
	    
    double getTargetIonDensityProfile(double,double,double);
    double getElectronPressureProfile(double,double,double);
    
    double getDensity(double,double,double,int);
    std::vector<double> getVelocity(double,double,double,int);
    std::vector<double> getBfield(double,double,double);
   
    std::vector<double> getFluidVelocity(double,double,double,int);
    std::vector<double> getFluidVelocity4InjectedParticles(double,double,double);
    std::vector<double> getVelocity4InjectedParticles(double,double,double);
    double getTimeStep();
    int getMaxTimestepsNum();
    int getNumberOfSpecies();
    int getTimestepsNum2Write();

    double getCollisionFrequencyFactor();
    double getDefaultCoulombLogarithm();
    
    double getPPC4species(int);
    double getMass4species(int);
    double getCharge4species(int);
    int getIfSpeciesFrozen(int);

    int getDFtype(int);
    int getDFtype4InjectedParticles();
    
    std::string getOutputDir() const;
    std::string getFilenameTemplate();
    
    PyObject * getPythonClassInstance(std::string className);
};
#endif /* Loader_hpp */
