#ifndef EleMagManager_hpp
#define EleMagManager_hpp

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


class EleMagManager{
    
private:
    std::unique_ptr<Logger> logger;
    std::shared_ptr<Loader> loader;
    std::shared_ptr<GridManager> gridMgr;
    
    double* BfieldDampingCoeff;
    void initBfieldDampingCoeff();
    
    void initialize();
    void calculateMagneticField(int, int, int);
    void calculateCurrent(int, int);
    
    void write2Log(int, int, int, int, const double*,
                   double*, double*, double*, const double*, double);
    
public:
    EleMagManager(std::shared_ptr<Loader>, std::shared_ptr<GridManager>);
    ~EleMagManager();
    void calculateBhalf(int);
    void calculateBnext();
    void calculateJhalf(int);
    void calculateJnext();
    
    void calculateEnext(int);

};


#endif /* EleMagManager_hpp */
