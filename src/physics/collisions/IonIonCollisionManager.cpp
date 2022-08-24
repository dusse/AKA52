#include "IonIonCollisionManager.hpp"

using namespace std;
using namespace chrono;

IonIonCollisionManager::IonIonCollisionManager(shared_ptr<Loader> ldr,
                                               shared_ptr<GridManager> gridMnr,
                                               shared_ptr<Pusher> pshr):
                                               loader(move(ldr)), 
                                               gridMgr(move(gridMnr)),
                                               pusher(move(pshr))
{


    logger.reset(new Logger());
    initialize();
    logger->writeMsg("[IonIonCollisionManager] create...OK", DEBUG);
}

void IonIonCollisionManager::initialize(){
}

IonIonCollisionManager::~IonIonCollisionManager(){
}

void IonIonCollisionManager::collideIons(int phase){
    auto start_time = high_resolution_clock::now();
    string msg0 ="[IonIonCollisionManager] start to collide ions ";
    logger->writeMsg(msg0.c_str(), DEBUG);

    int velShift = 0;
    if( phase == CORRECTOR ){
        velShift = 3;
    }

    int numOfSpecies = loader->getNumberOfSpecies();
    Particle** particles = pusher->getParticles();
    int totalPrtclNumber = pusher->getTotalParticleNumber();

    int type, type2, idx, idxG2;
    
    int xSize = loader->resolution[0];
    int ySize = loader->resolution[1];
    int zSize = loader->resolution[2];
    
    int xSizeG2 = xSize+2;
    int ySizeG2 = ySize+2;
    int zSizeG2 = zSize+2;
    
    int G2nodesNumber = xSizeG2*ySizeG2*zSizeG2;

    int* particlesNumber = new int[G2nodesNumber*numOfSpecies];
    map<int, map<int, vector<int>>> particlesInEachCell;

    for( idx = 0; idx < G2nodesNumber; idx++ ){
        map<int, vector<int>> particlesOfTheGivenType;
        vector<int> particleIndecies;
        for( type = 0; type < numOfSpecies; type++ ){
            particlesNumber[numOfSpecies*idx+type] = 0;
            particlesOfTheGivenType[type] = particleIndecies;
        }
        particlesInEachCell[idx] = particlesOfTheGivenType;
    }
    
    double* pos;
    double* vel;
    double pw, mass;
    map<int, VectorVar**> dens_vel;
    for( type = 0; type < numOfSpecies; type++ ){
            dens_vel[type] = gridMgr->getVectorVariableOnG2(gridMgr->DENS_VEL(type));
    }

    int ptclIdx, ion1idx, ion2idx;
    double dens1, dens2;
    double pw1, pw2;
    int group1Idx, group2Idx, restIdx;
    for( idxG2 = 0; idxG2 < G2nodesNumber; idxG2++ ){

        /*** shuffle particles ***/
        for( type = 0; type < numOfSpecies; type++ ){
            random_shuffle(particlesInEachCell[idxG2][type].begin(), particlesInEachCell[idxG2][type].end());
        }
        /** based on work of Nicolas Loic (2017)
         * Effects of collisions on the magnetic streaming instability 
         * (Doctoral dissertation, Paris 6).**/        
        for( type = 0; type < numOfSpecies; type++ ){
            if( pusher->getIfParticleTypeIsFrozen(type) == 1 ) continue;

            dens1 = dens_vel[type][idxG2]->getValue()[0];

            int numOfPartclsOfGvnType = particlesNumber[numOfSpecies*idxG2+type];
            if( numOfPartclsOfGvnType == 0 ) continue;
            /** intra-species collisions
                Takizuka, T. and H. Abe (1977). “A binary collision model for plasma simulation with a
                particle code” Journal of Computational Physics 25 **/
            if( numOfPartclsOfGvnType%2 == 0 ){
                for( ptclIdx = 0; ptclIdx < numOfPartclsOfGvnType/2; ptclIdx++ ){
                    ion1idx = particlesInEachCell[idxG2][type][2*ptclIdx];  
                    ion2idx = particlesInEachCell[idxG2][type][2*ptclIdx+1];
                    scatterVelocities(velShift, ion1idx, ion2idx, type, type, dens1, dens1, dens1,
                                      particles[ion1idx]->getVelocity(), particles[ion2idx]->getVelocity(),
                                      1);
                }
            }else{
                for( ptclIdx = 0; ptclIdx < (numOfPartclsOfGvnType/2)-1; ptclIdx++ ){
                    ion1idx = particlesInEachCell[idxG2][type][2*ptclIdx];
                    ion2idx = particlesInEachCell[idxG2][type][2*ptclIdx+1];
                    scatterVelocities(velShift, ion1idx, ion2idx, type, type, dens1, dens1, dens1,
                                      particles[ion1idx]->getVelocity(), particles[ion2idx]->getVelocity(),
                                      1);
                }
                if( numOfPartclsOfGvnType >= 3 ){
                    ion1idx = particlesInEachCell[idxG2][type][numOfPartclsOfGvnType-2];
                    ion2idx = particlesInEachCell[idxG2][type][numOfPartclsOfGvnType-1];
                    scatterVelocities(velShift, ion1idx, ion2idx, type, type, dens1, dens1, dens1,
                                      particles[ion1idx]->getVelocity(), particles[ion2idx]->getVelocity(),
                                      0.5);

                    ion1idx = particlesInEachCell[idxG2][type][numOfPartclsOfGvnType-3];
                    ion2idx = particlesInEachCell[idxG2][type][numOfPartclsOfGvnType-1];
                    scatterVelocities(velShift, ion1idx, ion2idx, type, type, dens1, dens1, dens1,
                                      particles[ion1idx]->getVelocity(), particles[ion2idx]->getVelocity(),
                                      0.5);

                    ion1idx = particlesInEachCell[idxG2][type][numOfPartclsOfGvnType-3];
                    ion2idx = particlesInEachCell[idxG2][type][numOfPartclsOfGvnType-2];
                    scatterVelocities(velShift, ion1idx, ion2idx, type, type, dens1, dens1, dens1,
                                      particles[ion1idx]->getVelocity(), particles[ion2idx]->getVelocity(),
                                      0.5);
                }
            }
            
            pw1 = pusher->getParticleWeight4Type(type);
            /** inter-species collisions
                Miller, R. H. and M. R. Combi (1994). “A Coulomb collision algorithm for weighted particle
                simulations” Geophysical Research Letters 21 **/
            for( type2 = type+1; type2 < numOfSpecies; type2++ ){
                if( pusher->getIfParticleTypeIsFrozen(type2) == 1 ) continue;

                int numOfPartclsOfGvnType2 = particlesNumber[numOfSpecies*idxG2+type2];
                if( numOfPartclsOfGvnType2 == 0 ) continue;

                pw2 = pusher->getParticleWeight4Type(type2);
                dens2 = dens_vel[type2][idxG2]->getValue()[0];
                
                if( numOfPartclsOfGvnType >= numOfPartclsOfGvnType2 ){

                    double weightFactor = (pw1 >= pw2) ? pw1/pw2 : 1;

                    int    quotient  = floor((double)numOfPartclsOfGvnType/(double)numOfPartclsOfGvnType2);
                    double remainder = ((double)numOfPartclsOfGvnType/(double)numOfPartclsOfGvnType2)-quotient;

                    int firstGroupSpecie2 = (int)round(remainder*double(numOfPartclsOfGvnType2));

                    for( group1Idx = 0; group1Idx < firstGroupSpecie2; group1Idx++ ){
                        ion1idx = particlesInEachCell[idxG2][type2][group1Idx];
                        for (restIdx = 0; restIdx < quotient+1; restIdx++) {
                            ptclIdx = group1Idx*(quotient+1)+restIdx;
                            ion2idx = particlesInEachCell[idxG2][type][ptclIdx];
                            scatterVelocities(velShift, ion1idx, ion2idx, type2, type, dens2, dens1,
                                         min(dens1,dens2), particles[ion1idx]->getVelocity(), particles[ion2idx]->getVelocity(),
                                         weightFactor);
                        }
                    }

                    for( group2Idx = firstGroupSpecie2; group2Idx < numOfPartclsOfGvnType2; group2Idx++ ){
                        ion1idx = particlesInEachCell[idxG2][type2][group2Idx];
                        for( restIdx = 0; restIdx < quotient; restIdx++ ){
                            ptclIdx = firstGroupSpecie2*(quotient+1)+(group2Idx-firstGroupSpecie2)*quotient+restIdx;
                            ion2idx = particlesInEachCell[idxG2][type][ptclIdx];
                            scatterVelocities(velShift, ion1idx, ion2idx, type2, type, dens2, dens1,
                                         min(dens1,dens2), particles[ion1idx]->getVelocity(), particles[ion2idx]->getVelocity(),
                                         weightFactor);
                        }
                    }
                }else{
                    double weightFactor = (pw1 >= pw2) ? 1 : pw2/pw1;

                    int    quotient  = floor((double)numOfPartclsOfGvnType2/(double)numOfPartclsOfGvnType);
                    double remainder = ((double)numOfPartclsOfGvnType2/(double)numOfPartclsOfGvnType)-quotient;
                    int firstGroupSpecie1 = (int)round(remainder*numOfPartclsOfGvnType);

                    for( group1Idx = 0; group1Idx < firstGroupSpecie1; group1Idx++ ){
                        ion1idx = particlesInEachCell[idxG2][type][group1Idx];
                        for( restIdx = 0; restIdx < quotient+1; restIdx++ ){
                            ptclIdx = group1Idx*(quotient+1)+restIdx;
                            ion2idx = particlesInEachCell[idxG2][type2][ptclIdx];
                            scatterVelocities(velShift, ion1idx, ion2idx, type, type2, dens1, dens2,
                                         min(dens1,dens2), particles[ion1idx]->getVelocity(), particles[ion2idx]->getVelocity(),
                                         weightFactor);
                        }
                    }

                    for( group2Idx = firstGroupSpecie1; group2Idx < numOfPartclsOfGvnType; group2Idx++ ){
                        ion1idx = particlesInEachCell[idxG2][type][group2Idx];
                        for( restIdx = 0; restIdx < quotient; restIdx++ ){
                            ptclIdx = firstGroupSpecie1*(quotient+1)+(group2Idx-firstGroupSpecie1)*quotient+restIdx;
                            ion2idx = particlesInEachCell[idxG2][type2][ptclIdx];
                            scatterVelocities(velShift, ion1idx, ion2idx, type, type2, dens1, dens2,
                                         min(dens1,dens2), particles[ion1idx]->getVelocity(), particles[ion2idx]->getVelocity(),
                                         weightFactor);
                        }
                    }
                }
            }
        }
    }


    delete [] particlesNumber;

    auto end_time = high_resolution_clock::now();
    auto msg ="[IonIonCollisionManager] collideIons()... duration = "
                +to_string(duration_cast<milliseconds>(end_time - start_time).count())+" ms";
    logger->writeMsg(msg.c_str(), DEBUG);
}

void IonIonCollisionManager::scatterVelocities(int velShift, int ion1idx, int ion2idx, int type1, int type2,
                                               double dens1, double dens2, double lowestDensity,
                                               double* ion1Vel, double* ion2Vel, double factor ){

    double charge1 = pusher->getParticleCharge4Type(type1);
    double charge2 = pusher->getParticleCharge4Type(type2);
    double mass1 = pusher->getParticleMass4Type(type1);
    double mass2 = pusher->getParticleMass4Type(type2);
    double reducedMass = mass1*mass2/(mass1+mass2);
    double prtcleWeight1 = pusher->getParticleWeight4Type(type1);
    double prtcleWeight2 = pusher->getParticleWeight4Type(type2);
    double relativeVelX = ion1Vel[0+velShift]-ion2Vel[0+velShift];
    double relativeVelY = ion1Vel[1+velShift]-ion2Vel[1+velShift];
    double relativeVelZ = ion1Vel[2+velShift]-ion2Vel[2+velShift];
    double relativeVelPerp = sqrt(pow(relativeVelX,2)+pow(relativeVelY,2));
    double relativeVelMod  = sqrt(pow(relativeVelX,2)+pow(relativeVelY,2)+pow(relativeVelZ,2));
    relativeVelMod = (relativeVelMod > EPSILON) ? relativeVelMod : EPSILON;
    double collisionFrequency;
    double rndm1, rndm2, rndm3;
    double defaultCollisionFrequencyFactor = factor*loader->getCollisionFrequencyFactor()*loader->getTimeStep();
    double coulombLog = loader->getDefaultCoulombLogarithm();

    collisionFrequency = (defaultCollisionFrequencyFactor*pow(charge1,2)*pow(charge2,2)*lowestDensity*coulombLog)
                            /(pow(reducedMass,2)*pow(relativeVelMod,3));

    rndm1 = RNM;
    rndm1 = -2*log((rndm1 > EPSILON) ? rndm1 : EPSILON);
    rndm2 = 2*PI*RNM;

    double variance = sqrt(collisionFrequency*rndm1)*cos(rndm2);
    double sinTheta, oneMinusCosTheta, phi;

    if( collisionFrequency < 1 ){
        sinTheta = 2*variance/(1+pow(variance,2));
        oneMinusCosTheta = 2*pow(variance,2)/(1+pow(variance,2));
    }else{
        rndm3 = 2*PI*RNM;
        sinTheta = sin(rndm3);
        oneMinusCosTheta = 1-cos(rndm3);
    }
    phi = 2*PI*RNM;

    double alpha, betta;
    rndm1 = RNM;
    alpha = (rndm1 < prtcleWeight2/prtcleWeight1) ? 1.0 : 0.0;
    betta = (rndm1 < prtcleWeight1/prtcleWeight2) ? 1.0 : 0.0;
    double velDeltas[3] = {0.0,0.0,0.0};
    if( relativeVelPerp == 0.0 ){
        velDeltas[0] =  relativeVelMod*sinTheta*cos(phi);
        velDeltas[1] =  relativeVelMod*sinTheta*sin(phi);
        velDeltas[2] = -relativeVelMod*oneMinusCosTheta;
    }else{
        double velXnorm = relativeVelX/relativeVelPerp;
        double velYnorm = relativeVelY/relativeVelPerp;
        velDeltas[0] = velXnorm*relativeVelZ*sinTheta*cos(phi)
                      -velYnorm*relativeVelMod*sinTheta*sin(phi)
                      -relativeVelX*oneMinusCosTheta;

        velDeltas[1] = velYnorm*relativeVelZ*sinTheta*cos(phi)
                      +velXnorm*relativeVelMod*sinTheta*sin(phi)
                      -relativeVelY*oneMinusCosTheta;

        velDeltas[2] = -relativeVelPerp*sinTheta*cos(phi)
                       -relativeVelZ*oneMinusCosTheta;
    }

    for( int coord = 0; coord < 3; coord++ ){
        pusher->setParticleVelocity(ion1idx, coord+velShift, 
                                    ion1Vel[coord+velShift]+alpha*(reducedMass/mass1)*velDeltas[coord]);
        pusher->setParticleVelocity(ion2idx, coord+velShift,
                                    ion2Vel[coord+velShift]-betta*(reducedMass/mass2)*velDeltas[coord]);
    }
}
