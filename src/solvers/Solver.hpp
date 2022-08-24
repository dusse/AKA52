//
//  Solver.hpp

#ifndef Solver_hpp
#define Solver_hpp

#include <chrono>
#include <stdio.h>
#include <iostream>
#include <string>
#include <memory>

#include "../grid/GridManager.hpp"
#include "../input/Loader.hpp"
#include "../misc/Misc.hpp"

#include "../physics/pusher/Pusher.hpp"
#include "../physics/hydro/HydroManager.hpp"
#include "../physics/electro-magnetic/EleMagManager.hpp"
#include "../physics/pressure-closure/ClosureManager.hpp"
#include "../physics/laser/LaserMockManager.hpp"
#include "../physics/collisions/IonIonCollisionManager.hpp"


const static int  SOLVE_OK   = 0;
const static int  SOLVE_FAIL = 1;

class Solver
{
    
private:
    
    std::unique_ptr<Logger> logger;
    
    std::shared_ptr<Loader> loader;
    std::shared_ptr<GridManager> gridMng;
    std::shared_ptr<Pusher> pusher;
    std::shared_ptr<HydroManager> hydroMng;
    std::shared_ptr<EleMagManager> emMng;
    std::shared_ptr<ClosureManager> closureMng;
    std::shared_ptr<LaserMockManager> laserMng;
    std::shared_ptr<IonIonCollisionManager> collideMng;
    
    
    void performCalculation(int, int);

    
public:
    Solver(std::shared_ptr<Loader>, std::shared_ptr<GridManager>,
                 std::shared_ptr<Pusher>, std::shared_ptr<HydroManager>,
                 std::shared_ptr<EleMagManager>, std::shared_ptr<ClosureManager>,
                 std::shared_ptr<LaserMockManager>, std::shared_ptr<IonIonCollisionManager>);
        
    void initialize();
    int solve(int);
    void finilize();
    ~Solver();
};
#endif
