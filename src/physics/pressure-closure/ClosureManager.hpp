#ifndef ClosureManager_hpp
#define ClosureManager_hpp

#include <stdio.h>

#include <chrono>
#include <stdio.h>
#include <iostream>
#include <cmath>
#include <string>
#include <memory>

#include "../../grid/GridManager.hpp"
#include "../../input/Loader.hpp"
#include "../../misc/Misc.hpp"

class ClosureManager{
    
private:
    std::unique_ptr<Logger> logger;
    std::shared_ptr<Loader> loader;
    std::shared_ptr<GridManager> gridMgr;
    
    double emass = 0.01;
    
    double* driverNext;
    double* pressuNext;
    double* electrnVel;
    double* bfieldPrev;
    double* pressuInit;
    double* pressureDampingCoeff;
    
    void initialize();
    void initPressure();
    void initPressureDampingCoeff();
    
    void subCycledPressure(int, int);
    void implicitPressure(int, int);
    
    void setDriver(int );
    void setIsotropization(double[6], double[6]);
    
    void gradients(double[3][3], double[3][3], double[3][3][3], int[3]);
    
    void transformMatrix(double[3][3], double[3][3], double[3][3], int);
    void ortho(double[3], double[3][3]);
    
public:
    ClosureManager(std::shared_ptr<Loader>, std::shared_ptr<GridManager>);
    ~ClosureManager();
    void calculatePressure(int, int);

};


#endif /* ClosureManager_hpp */
