#ifndef HydroManager_hpp
#define HydroManager_hpp

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


class HydroManager{
    
private:
    std::unique_ptr<Logger> logger;
    std::shared_ptr<Loader> loader;
    std::shared_ptr<GridManager> gridMgr;
    std::shared_ptr<Pusher> pusher;
    std::shared_ptr<BoundaryManager> boundaryMgr;
    
    void initialize();
    void calculateAvgFluidVelocity4AllSpiecies(int);
    
public:
    HydroManager(std::shared_ptr<Loader>, std::shared_ptr<GridManager>, std::shared_ptr<Pusher>);
    void gatherMoments(int);
    

};
#endif /* HydroManager_hpp */
