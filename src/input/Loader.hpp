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
    
    double hyperviscosity;

    double minimumDens2ResolvePPC;
    
    int smoothStride;
    double electronmass;
    double relaxFactor;
    
    //stability parameters
    double cellBreakdownEfieldFactor = 0.005;
    double criticalPressure = 500.0;
    
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
    
    std::vector<double> getVelocity(double,double,double,int);
    std::vector<double> getBfield(double,double,double);
    double getDensity(double,double,double,int);
    double getTimeStep();
    int getMaxTimestepsNum();
    int getNumberOfSpecies();
    int getTimestepsNum2Write();
    double getMass4spieceies(int);
    double getCharge4spieceies(int);
    double getParticlesPerCellNumber();
    std::string getOutputDir() const;
    std::string getFilenameTemplate();
    PyObject * getPythonClassInstance(std::string className);
};
#endif /* Loader_hpp */
