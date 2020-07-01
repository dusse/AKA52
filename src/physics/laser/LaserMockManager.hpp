#ifndef LaserMockManager_hpp
#define LaserMockManager_hpp

#include <stdio.h>

#include <chrono>
#include <stdio.h>
#include <iostream>
#include <cmath>
#include <string>
#include <memory>
#include<algorithm>

#include "../../grid/GridManager.hpp"
#include "../../particles/Particle.hpp"
#include "../pusher/Pusher.hpp"
#include "../../input/Loader.hpp"
#include "../../misc/Misc.hpp"
#include "../../common/variables/VectorVar.hpp"

class LaserMockManager{
    
private:
    std::unique_ptr<Logger> logger;
    std::shared_ptr<Loader> loader;
    std::shared_ptr<GridManager> gridMgr;
    std::shared_ptr<Pusher> pusher;
    
    int PARTICLE_TYPE_HEAT = 1;
    double PARTICLE_DENS2KEEP = 2.0;
    double PARTICLE_TEMP2KEEP = 1.0;
    
    void initialize();

    int partTemp(double[3],  double[2]);
    int eleTemp(double[3],  double[2]);
    double calcDens(double[3]);
    
public:
    LaserMockManager(std::shared_ptr<Loader>,
                     std::shared_ptr<GridManager>,
                     std::shared_ptr<Pusher>);
    void addIons();
    void accelerate();
    

};


#endif /* LaserMockManager_hpp */
