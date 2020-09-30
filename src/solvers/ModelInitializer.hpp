//
//  ModelInitializer.hpp

#ifndef ModelInitializer_hpp
#define ModelInitializer_hpp

#include <chrono>
#include <stdio.h>
#include <iostream>
#include <string>
#include <memory>
#include <hdf5.h>

#include "../grid/GridManager.hpp"
#include "../input/Loader.hpp"
#include "../misc/Misc.hpp"
#include "../physics/pusher/Pusher.hpp"


class ModelInitializer
{
    
private:
    
    std::unique_ptr<Logger> logger;
    
    std::shared_ptr<Loader> loader;
    std::shared_ptr<GridManager> gridMng;
    std::shared_ptr<Pusher> pusher;
    
    double totalDensityInTHEbox;
    double globalMinimumDens;
    int totParticleNumberInDomain;
    
    std::vector<int> npc;
    std::vector<double> prtcleWeight;
    
    void initMagneticField();
    void initVariablesonG2();
    void initParticles();
    
    void readAllFromFile();
    void readField(hid_t, std::string, hid_t, hid_t, hid_t, void * );

    
public:
    ModelInitializer(std::shared_ptr<Loader>,
                     std::shared_ptr<GridManager>,
                     std::shared_ptr<Pusher>);
        
    void initialize();
    
    void finilize();
    ~ModelInitializer();
};
#endif
