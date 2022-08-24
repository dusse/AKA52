#ifndef IonIonCollisionManager_hpp
#define IonIonCollisionManager_hpp

#include <stdio.h>

#include <chrono>
#include <stdio.h>
#include <iostream>
#include <cmath>
#include <string>
#include <memory>

#include "../../grid/GridManager.hpp"
#include "../pusher/Pusher.hpp"
#include "../../input/Loader.hpp"
#include "../../misc/Misc.hpp"


class IonIonCollisionManager{
    
private:
    std::unique_ptr<Logger> logger;
    std::shared_ptr<Loader> loader;
    std::shared_ptr<GridManager> gridMgr;
    std::shared_ptr<Pusher> pusher;
    
    
    void initialize();
    
public:
    IonIonCollisionManager(std::shared_ptr<Loader>,
                           std::shared_ptr<GridManager>,
                           std::shared_ptr<Pusher>);
    ~IonIonCollisionManager();
    void collideIons(int);

    void scatterVelocities(int, int, int, int, int,
                           double, double ,double,
                           double*, double* , double);
};


#endif /* IonIonCollisionManager_hpp */
