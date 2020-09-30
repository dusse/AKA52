#ifndef Particle_hpp
#define Particle_hpp

#include <memory>
#include <stdio.h>
#include <iostream>
#include <string>
#include <map>

// number of double fields for MPI communication
const int PARTICLES_SIZE = 13;

class Particle{
    
private:
    // Predictor [0,1,2]  Corrector [3,4,5]
    double pos[6] = {0,0,0,0,0,0};
    double vel[6] = {0,0,0,0,0,0};
    
    int type = 0;

public:
    Particle();
    Particle(std::shared_ptr<Particle>);
    
    double* getPosition();
    void    setPosition(double[6]);
    void    setPosition(int, double);
    
    double* getVelocity();
    void    setVelocity(double[6]);
    void    setVelocity(int, double);
    
    int  getType();
    void setType(int);
    
    void serialize(double*, int);
    void deserialize(double *, int);
    
    void reinitializeUsingParticle(std::shared_ptr<Particle>);
    void reinitializeUsingParticle(Particle*);
    
};
#endif