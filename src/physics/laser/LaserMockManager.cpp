 #include "LaserMockManager.hpp"

using namespace std;
using namespace chrono;


LaserMockManager::LaserMockManager(std::shared_ptr<Loader> ldr,
                           std::shared_ptr<GridManager> gridMnr,
                           shared_ptr<Pusher> pshr):loader(move(ldr)),
                            gridMgr(move(gridMnr)), pusher(move(pshr)){
                                
    logger.reset(new Logger());
    initialize();
    logger->writeMsg("[LaserMockManager] create...OK", DEBUG);
}


LaserMockManager::~LaserMockManager(){
    delete[] electronPressureProfile;
    delete[] targetIonDensityProfile;
    delete[] ionThermalVelocityProfile;
    delete[] ionFluidVelocityProfile;
}


void LaserMockManager::initialize(){
    int numOfSpecies = loader->getNumberOfSpecies();

    PARTICLE_TYPE2LOAD = numOfSpecies-1;//last type is reserved for loaded particles
    
    int xRes = loader->resolution[0],
        yRes = loader->resolution[1],
        zRes = loader->resolution[2];
    
    int xResG2 = xRes+2, yResG2 = yRes+2, zResG2 = zRes+2;
    
    int nG2 = xResG2*yResG2*zResG2;
    
    electronPressureProfile = new double[nG2*sizeof(double)];
    targetIonDensityProfile = new double[nG2*sizeof(double)];
    ionThermalVelocityProfile = new double[nG2*3*sizeof(double)];
    ionFluidVelocityProfile = new double[nG2*3*sizeof(double)];

    double pres, dens;
    vector<double> fluidVel;
    vector<double> thermalVel;

    int idx, idxOnG2;
    
    double domainShiftX = loader->boxCoordinates[0][0];
    double domainShiftY = loader->boxCoordinates[1][0];
    double domainShiftZ = loader->boxCoordinates[2][0];
    
    double dx = loader->spatialSteps[0];
    double dy = loader->spatialSteps[1];
    double dz = loader->spatialSteps[2];
    
    double energyOnDomain = 0.0;
    double G2shift = -0.5;// symmetry requirement
    double x,y,z;
    int i,j,k,n;
    for( i = 0; i < xRes+1; i++){
        for( j = 0; j < yRes+1; j++) {
            for( k = 0; k < zRes+1; k++){
                
                idxOnG2 = IDX(i, j, k, xResG2, yResG2, zResG2);
                
                x = (i + G2shift)*dx + domainShiftX;
                y = (j + G2shift)*dy + domainShiftY;
                z = (k + G2shift)*dz + domainShiftZ;
                
                pres = loader->getElectronPressureProfile(x,y,z);
                
                electronPressureProfile[idxOnG2] = pres;
                
                energyOnDomain += pres;
                
                dens = loader->getTargetIonDensityProfile(x,y,z);
                
                targetIonDensityProfile[idxOnG2] = dens;
                
                thermalVel = loader->getVelocity4InjectedParticles(x, y, z);
                fluidVel = loader->getFluidVelocity4InjectedParticles(x, y, z);

                for( n = 0; n < 3; n++ ){
                    ionThermalVelocityProfile[3*idxOnG2+n] = thermalVel[n];
                    ionFluidVelocityProfile[3*idxOnG2+n] = fluidVel[n];
                }                
            }
        }
    }
    
    double TOT_IN_BOX = 0;
    MPI_Allreduce(&energyOnDomain, &TOT_IN_BOX, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    logger->writeMsg(("[LaserMockManager] total energy loaded in the box per step = "
                      +to_string(TOT_IN_BOX*loader->getTimeStep())).c_str(),  DEBUG);
    
}


void LaserMockManager::addIons(int i_time){
    auto start_time = high_resolution_clock::now();

    const double VELOCITY4COLDTEMPERATURE = EPSILON;

    double dx = loader->spatialSteps[0];
    double dy = loader->spatialSteps[1];
    double dz = loader->spatialSteps[2];
    double dl[3] = {dx, dy, dz};
    
    double domainShiftX = loader->boxCoordinates[0][0];
    double domainShiftY = loader->boxCoordinates[1][0];
    double domainShiftZ = loader->boxCoordinates[2][0];
    
    int xRes = loader->resolution[0],
        yRes = loader->resolution[1],
        zRes = loader->resolution[2];
    int xResG2 = xRes+2, yResG2 = yRes+2, zResG2 = zRes+2;
        
    int i,j,k;
    double vpb[3] = {0.0, 0.0, 0.0};
    double pos[3] = {0.0, 0.0, 0.0};
    
    int ptclIDX;
    
    VectorVar** dens4Injected = gridMgr->getVectorVariableOnG2(gridMgr->DENS_VEL(PARTICLE_TYPE2LOAD));
    VectorVar** dens4nonInjected = gridMgr->getVectorVariableOnG2(gridMgr->DENS_VEL(loader->prtclType2Load));

    vector<shared_ptr<Particle>> particles2add;
    int particle_idx = 0;
    double r1, r2;
    int idxOnG2;
    int requiredPrtclNum;
    int distributionType = loader->getDFtype4InjectedParticles();
    for( i = 0; i < xRes; i++ ){
        for( j = 0; j < yRes; j++ ){
            for( k = 0; k < zRes; k++ ){
                
                idxOnG2 = IDX(i+1, j+1, k+1, xResG2, yResG2, zResG2);
                
                double desireDens = 0.0;
                
                double pres = electronPressureProfile[idxOnG2];
                int type2use = 0;
                if( pres > 0.0 ){
                    type2use = PARTICLE_TYPE2LOAD;
                    desireDens = targetIonDensityProfile[idxOnG2] - dens4Injected[idxOnG2]->getValue()[0];
                }else{
                    type2use = loader->prtclType2Load;
                    desireDens = targetIonDensityProfile[idxOnG2] - dens4nonInjected[idxOnG2]->getValue()[0];
                }
                double particleWeight = pusher->getParticleWeight4Type(type2use);
                
                requiredPrtclNum = int(desireDens/particleWeight);
                
                for( ptclIDX = 0; ptclIDX < requiredPrtclNum; ptclIDX++ ){
                    
                    particles2add.push_back(shared_ptr<Particle>(new Particle));
                    
                    pos[0] = (i + RNM) * dx;
                    pos[1] = (j + RNM) * dy;
                    pos[2] = (k + RNM) * dz;
                    
                    pos[0] = (pos[0] == xRes*dx) ? xRes*dx - EPS4 : pos[0];
                    pos[1] = (pos[1] == yRes*dy) ? yRes*dy - EPS4 : pos[1];
                    pos[2] = (pos[2] == zRes*dz) ? zRes*dz - EPS4 : pos[2];
                    
                    pos[0] += domainShiftX;
                    pos[1] += domainShiftY;
                    pos[2] += domainShiftZ;
                    
                    double pos2Save[6] = {pos[0], pos[1], pos[2],
                                          pos[0], pos[1], pos[2]};
                    
                    particles2add[particle_idx]->setPosition(pos2Save);
                    particles2add[particle_idx]->setType(type2use);
                        
                    if( distributionType == 0 ){
                        r1 = RNM;
                        r2 = RNM;
                        r1   = (fabs(r1 - 1.0) < EPS8) ? r1 - EPS8 : r1;
                        r1   = (r1 > EPS8)? r1 : r1 + EPS8;
                        if( i_time > loader->laserPulseDuration_tsnum ){
                            vpb[0] = sqrt(-2*log(r1))*VELOCITY4COLDTEMPERATURE*cos(2*PI*r2);
                            vpb[1] = sqrt(-2*log(r1))*VELOCITY4COLDTEMPERATURE*sin(2*PI*r2);
                            r1 = RNM;
                            r2 = RNM;
                            r1   = (fabs(r1 - 1.0) < EPS8) ? r1 - EPS8 : r1;
                            r1   = (r1 > EPS8)? r1 : r1 + EPS8;
                            vpb[2] = sqrt(-2*log(r1))*VELOCITY4COLDTEMPERATURE*cos(2*PI*r2);
                        }else{
                            vpb[0] = sqrt(-2*log(r1))*ionThermalVelocityProfile[3*idxOnG2+0]*cos(2*PI*r2);
                            vpb[0]+= ionFluidVelocityProfile[3*idxOnG2+0];
                            vpb[1] = sqrt(-2*log(r1))*ionThermalVelocityProfile[3*idxOnG2+1]*sin(2*PI*r2);
                            vpb[1]+= ionFluidVelocityProfile[3*idxOnG2+1];
                            r1 = RNM;
                            r2 = RNM;
                            r1   = (fabs(r1 - 1.0) < EPS8) ? r1 - EPS8 : r1;
                            r1   = (r1 > EPS8)? r1 : r1 + EPS8;
                            vpb[2] = sqrt(-2*log(r1))*ionThermalVelocityProfile[3*idxOnG2+2]*cos(2*PI*r2);
                            vpb[2]+= ionFluidVelocityProfile[3*idxOnG2+2];
                         }

                    }else{
                        double randoms[3] = {RNM,RNM,RNM};
                        for( int comp = 0; comp < 3; comp++ ){
                            if( i_time > loader->laserPulseDuration_tsnum ){
                                vpb[comp] = (1.0-2.0*randoms[comp])*VELOCITY4COLDTEMPERATURE;
                            }else{
                                vpb[comp] = (1.0-2.0*randoms[comp])*ionThermalVelocityProfile[3*idxOnG2+comp];
                                vpb[comp]+= ionFluidVelocityProfile[3*idxOnG2+comp];
                            }
                        }
                    }
                        
                    double vel2Save[6] = {vpb[0], vpb[1], vpb[2],
                                          vpb[0], vpb[1], vpb[2]};
                    particles2add[particle_idx]->setVelocity(vel2Save);
                   
                    particle_idx++;
                }
            }
        }
    }
    pusher->addParticles(particles2add);
    
    auto end_time = high_resolution_clock::now();
    string msg ="[LaserMockManager] addIons() added = "+to_string(particle_idx)
                +" particles, duration = "
                +to_string(duration_cast<milliseconds>(end_time - start_time).count())+" ms";
    logger->writeMsg(msg.c_str(), DEBUG);
}



void LaserMockManager::accelerate(int i_time){
    auto start_time = high_resolution_clock::now();
    
    if (i_time > loader->laserPulseDuration_tsnum){
        return;
    }
    int i, j, k, idxOnG2;
    
    int xRes = loader->resolution[0],
        yRes = loader->resolution[1],
        zRes = loader->resolution[2];    
    int xResG2 = xRes+2, yResG2 = yRes+2, zResG2 = zRes+2;
    int G2nodesNumber = xResG2*yResG2*zResG2;

    VectorVar** pdriverr    = gridMgr->getVectorVariableOnG2(DRIVER);
    double energyOnDomain = 0;
    for( idxOnG2 = 0; idxOnG2 < G2nodesNumber; idxOnG2++ ){
        const double* dr = pdriverr[idxOnG2]->getValue();
        energyOnDomain += (dr[0]+dr[3]+dr[5])/3;
    }
    
    for( i = 1; i < xRes + 1; i++ ){
        for( j = 1; j < yRes + 1; j++ ){
            for( k = 1; k < zRes + 1; k++ ){                
                idxOnG2 = IDX(i ,j ,k , xResG2, yResG2, zResG2);                
                double pres2set = electronPressureProfile[idxOnG2];
                gridMgr->addVectorVariableForNodeG2(idxOnG2, DRIVER, 0, pres2set);
                gridMgr->addVectorVariableForNodeG2(idxOnG2, DRIVER, 3, pres2set);
                gridMgr->addVectorVariableForNodeG2(idxOnG2, DRIVER, 5, pres2set);
            }
        }
    }

    gridMgr->sendBoundary2Neighbor(DRIVER);
    gridMgr->applyBC(DRIVER);
    
    double energyOnDomainNew = 0;
    for( idxOnG2 = 0; idxOnG2 < G2nodesNumber; idxOnG2++ ){
        const double* dr = pdriverr[idxOnG2]->getValue();
        energyOnDomainNew += (dr[0]+dr[3]+dr[5])/3;
    }

    double laserdelta = (energyOnDomainNew - energyOnDomain)*loader->getTimeStep();
    double TOT_IN_BOX = 0;
    MPI_Allreduce(&laserdelta, &TOT_IN_BOX, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    logger->writeMsg(("[LaserMockManager] instant added driver  = "
                      +to_string(laserdelta)).c_str(),  DEBUG);
    loader->loadedEnergyPerStep += laserdelta;


    auto end_time = high_resolution_clock::now();
    string msg ="[LaserMockManager] accelerate() duration = "
                +to_string(duration_cast<milliseconds>(end_time - start_time).count())+" ms";
    logger->writeMsg(msg.c_str(), DEBUG);
    
}



