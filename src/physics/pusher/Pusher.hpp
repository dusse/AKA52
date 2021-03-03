#ifndef Pusher_hpp
#define Pusher_hpp

#include <stdio.h>

#include <chrono>
#include <stdio.h>
#include <iostream>
#include <cmath>
#include <string>
#include <memory>

#include "../../grid/GridManager.hpp"
#include "../../grid/boundary/BoundaryManager.hpp"

#include "../../particles/Particle.hpp"

#include "../../input/Loader.hpp"
#include "../../misc/Misc.hpp"


class Pusher{
    
private:
    std::unique_ptr<Logger> logger;
    std::shared_ptr<Loader> loader;
    std::shared_ptr<GridManager> gridMgr;
    std::shared_ptr<BoundaryManager> boundaryMgr;
    
    double* weights;
    double* charges;
    double* masses;
    
    int totinBoxInit = 0;
    int totalNum;
    
    double INITIAL_B_FIELD = 0.0;
    
    int currentPartclNumOnDomain = 0;
    Particle** particles;
    
    void initialize();
    
    int checkParticle(int, Particle*, std::string);
    
    void performSorting();
    
    void reallocateParticles();
    
public:
    Pusher(std::shared_ptr<Loader>,
           std::shared_ptr<GridManager>,
           std::shared_ptr<BoundaryManager>);
    
    ~Pusher();
    void push(int, int);
    
    double getParticleWeight4Type(int);
    double getParticleCharge4Type(int);
    double getParticleMass4Type(int);
    int getTotalParticleNumber();
    void setTotalParticleNumber(int);
    
    void setParticleWeight4Type(int, double);
    void setParticleCharge4Type(int, double);
    void setParticleMass4Type(int, double);
    void initParticles(int,int);
    void addParticles(std::vector<std::shared_ptr<Particle>>);
    
    void setParticlePosition(int, double[6]);
    void setParticleVelocity(int, double[6]);
    void setParticleVelocity(int, int, double);
    void setParticleType(int, int);
    
    void checkEnergyBalance(int);
    
    Particle** getParticles();
    

};


#endif
