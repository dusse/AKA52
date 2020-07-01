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

void LaserMockManager::initialize(){
    PARTICLE_TYPE_HEAT = loader->prtclType2Heat;
    PARTICLE_DENS2KEEP = loader->prtclDens2Keep;
    PARTICLE_TEMP2KEEP = loader->prtclTemp2Keep;
}

void LaserMockManager::addIons(){
    auto start_time = high_resolution_clock::now();

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
    
    
    double mass = pusher->getParticleMass4Type(PARTICLE_TYPE_HEAT);
    double vth  = sqrt(loader->prtclTemp2Load/mass);
    
    int i,j,k;
    double vel[3] = {vth, vth, vth};
    double vpb[3] = {0.0, 0.0, 0.0};
    double pos[3] = {0.0, 0.0, 0.0};
    
    int ptclIDX;
    
    VectorVar** dens = gridMgr->getVectorVariableOnG2(gridMgr->DENS_VEL(PARTICLE_TYPE_HEAT));
    vector<shared_ptr<Particle>> particles2add;
    int particle_idx = 0;
    double r1, r2;
    int idxOnG2;
    int requiredPrtclNum;
    for( i = 0; i < xRes; i++){
        for( j = 0; j < yRes; j++) {
            for( k = 0; k < zRes; k++){
                
                idxOnG2 = IDX(i+1, j+1, k+1, xResG2, yResG2, zResG2);
                
                pos[0] = i * dx + domainShiftX;
                pos[1] = j * dy + domainShiftY;
                pos[2] = k * dz + domainShiftZ;
                
                double desireDens = calcDens(pos) - dens[idxOnG2]->getValue()[0];
                double particleWeight = pusher->getParticleWeight4Type(PARTICLE_TYPE_HEAT);
                requiredPrtclNum = int(desireDens/particleWeight);
                
                for (ptclIDX=0; ptclIDX < requiredPrtclNum; ptclIDX++){
                    
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
                        
                    particles2add[particle_idx]->setType(PARTICLE_TYPE_HEAT);
                        
                    r1 = RNM;
                    r2 = RNM;
                    r1   = (fabs(r1 - 1.0) < EPS8) ? r1 - EPS8 : r1;
                    r2   = (fabs(r2 - 1.0) < EPS8) ? r2 - EPS8 : r2;
                    r1   = (r1 > EPS8)? r1 : r1 + EPS8;
                    r2   = (r2 > EPS8)? r2 : r2 + EPS8;
                    vpb[0] = sqrt(-2*log(r1))*vel[0] * cos(2*PI*r2);
                    vpb[1] = sqrt(-2*log(r1))*vel[1] * sin(2*PI*r2);
                        
                    r1 = RNM;
                    r2 = RNM;
                    r1   = (fabs(r1 - 1.0) < EPS8) ? r1 - EPS8 : r1;
                    r2   = (fabs(r2 - 1.0) < EPS8) ? r2 - EPS8 : r2;
                    r1   = (r1 > EPS8)? r1 : r1 + EPS8;
                    r2   = (r2 > EPS8)? r2 : r2 + EPS8;
                    vpb[2] = sqrt(-2*log(r1))*vel[2] * cos(2*PI*r2);
                        
                        
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


int LaserMockManager::partTemp(double pos[3],  double T[2]){
    double T0, T1;
    double x, modX;
    double tx, ty, tz;
    double polynomVal[2] = {0.0,0.0};
    int flag = 0;
    
    int numOfSpots = loader->numOfSpots;
    double* centersOfSpots = loader->centersOfSpots;
    double* spotsRadius    = loader->spotsRadius;
    
    double A, B, C=0.2;
    int spotNum = 0;
    for (int h = 0; h < numOfSpots; h++){
        
        tx = pos[0] - centersOfSpots[3*h+0];
        ty = pos[1] - centersOfSpots[3*h+1];
        tz = pos[2] - centersOfSpots[3*h+2];
        
        x = 1.0;
        
        switch (loader->dim) {
            case 1:
                if(tx > 0.0){
                    continue;
                }
                A = spotsRadius[h];
                x = sqrt(pow(tx/A, 2));
                
                break;
            case 2:
                if(ty > 0.0){
                    continue;
                }
                A = spotsRadius[h];
                B = A;
                x = sqrt(pow(tx/A, 2)+pow(ty/C, 2));
                
                break;
            case 3:
                
                if(tz > 0.0){
                    continue;
                }
                A = spotsRadius[h];
                B = A;
                x = sqrt(pow(tx/A, 2)+pow(ty/B, 2)+pow(tz/C, 2));
                break;
                
            default:
                throw runtime_error("dimension is not 1/2/3");
        }
        
        
        if (x > 1.0) { x = 1.0; }
        
        if(x == 1.0 && flag == 0) {
            flag = 1;
            continue;
        }
        
        modX = fabs(x);
        
        polynomVal[h] =  -6.0*modX*modX*modX*modX*modX
        +15.0*x*x*x*x
        -10.0*modX*modX*modX
        +1;
        
        if(polynomVal[h] != 0.0){
            spotNum = h;
        }
    }
    
    T0 = PARTICLE_TEMP2KEEP*polynomVal[0];
    T1 = PARTICLE_TEMP2KEEP*polynomVal[1];
    T[0] = (T0+T1)/(2);
    T[1] = T[0];
    return spotNum;
}

int LaserMockManager::eleTemp(double pos[3],  double T[2]){
    double T0, T1;
    double x, modX;
    double tx, ty, tz;
    double polynomVal[2] = {0.0,0.0};
    int flag = 0;
    
    int numOfSpots = loader->numOfSpots;
    double* centersOfSpots = loader->centersOfSpots;
    double* spotsRadius    = loader->spotsRadius;
    
    double A, B, C=0.2;
    int spotNum = 0;
    for (int h = 0; h < numOfSpots; h++){
        
        tx = pos[0] - centersOfSpots[3*h+0];
        ty = pos[1] - centersOfSpots[3*h+1];
        tz = pos[2] - centersOfSpots[3*h+2];
        
        
        switch (loader->dim) {
            case 1:
                if(tx > 0.0){
                    continue;
                }
                A = spotsRadius[h];
                x = sqrt(pow(tx/A, 2));
                
                break;
            case 2:
                if(ty > 0.0){
                    continue;
                }
                A = spotsRadius[h];
                x = sqrt(pow(tx/A, 2)+pow(ty/C, 2));
                
                break;
            case 3:
                
                if(tz > 0.0){
                    continue;
                }
                A = spotsRadius[h];
                B = A;
                x = sqrt(pow(tx/A, 2)+pow(ty/B, 2)+pow(tz/C, 2));
                break;
                
            default:
                throw runtime_error("dimension is not 1/2/3");
        }
        
        if (x > 1.0) { x = 1.0; }
        
        if(x == 1.0 && flag == 0) {
            flag = 1;
            continue;
        }
        
        modX = fabs(x);
        
        polynomVal[h] =  -6.0*modX*modX*modX*modX*modX
        +15.0*x*x*x*x
        -10.0*modX*modX*modX
        +1;
        
        if(polynomVal[h] != 0.0){
            spotNum = h;
        }
    }
    
    T0 = PARTICLE_TEMP2KEEP*polynomVal[0];
    T1 = PARTICLE_TEMP2KEEP*polynomVal[1];
    T[0] = (T0+T1)/(2);
    T[1] = T[0];
    return spotNum;
}


double LaserMockManager::calcDens(double pos[3]){
    double x, modX;
    double tx, ty, tz;
    
    int flag = 0;
    int numOfSpots = loader->numOfSpots;
    double* centersOfSpots = loader->centersOfSpots;
    
    for (int h = 0; h < numOfSpots; h++){
        
        tx = pos[0] - centersOfSpots[3*h+0];
        ty = pos[1] - centersOfSpots[3*h+1];
        tz = pos[2] - centersOfSpots[3*h+2];
        switch (loader->dim) {
            case 1:
                if( tx > 0.0){
                    continue;
                }else{
                    return PARTICLE_DENS2KEEP;
                }
                
                break;
            case 2:
                if( ty > 0.0){
                    continue;
                }else{
                    return PARTICLE_DENS2KEEP;
                }
                
                break;
            case 3:
                
                if( tz > 0.0){
                    continue;
                }else{
                    return PARTICLE_DENS2KEEP;
                }
                
                break;
                
            default:
                throw runtime_error("dimension is not 1/2/3");
        }
        
    }
    return 0.0;
}


void LaserMockManager::accelerate(){
    auto start_time = high_resolution_clock::now();
    
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
    
    double xw, yw, zw;
    int i, j, k, l, m, idxOnG2;
    
    double pos[3], T[2], velRand[3], currentDens;
    const double factor = 1.0;
    
    double r1, r2;
    
    Particle** particles = pusher->getParticles();
    
    int totalPrtclNumber = pusher->getTotalParticleNumber();
    
    double * positio;
    double * velocit;
    int type;
    double mass;
    double predVel[3], corrVel[3];
    // update ions
    for (m = 0; m < totalPrtclNumber; m++){
        
        positio  = particles[m]->getPosition();
        velocit  = particles[m]->getVelocity();
        type = particles[m]->getType();
        
        if(type != PARTICLE_TYPE_HEAT){
            continue;
        }
        
        mass = pusher->getParticleMass4Type(type);
        
        pos[0] = positio[0];
        pos[1] = positio[1];
        pos[2] = positio[2];
        
        int spotNum = partTemp(pos, T);
        
        if (T[0] <= 0) {
            continue;
        };
        
        for(int coord = 0; coord<3; coord++){
            predVel[coord] = velocit[  coord];
            corrVel[coord] = velocit[3+coord];
        }
        
        double avgvel[3] = {0.0, 0.0, 0.0};
        
        avgvel[0] = 0.5*(predVel[0]+corrVel[0]);
        avgvel[1] = 0.5*(predVel[1]+corrVel[1]);
        avgvel[2] = 0.5*(predVel[2]+corrVel[2]);
        
        
        double temp = mass*(avgvel[0]*avgvel[0]
                           +avgvel[1]*avgvel[1]
                           +avgvel[2]*avgvel[2]);
        
        if(temp>T[0]){
            continue;
        }
        
        double V0 = sqrt(T[0]/mass);
        double VMAX = sqrt(PARTICLE_TEMP2KEEP/mass);
        double dV = abs(VMAX - V0);
        
        double* centersOfSpots = loader->centersOfSpots;
        
        double tx = pos[0] - centersOfSpots[3*spotNum+0];
        double ty = pos[1] - centersOfSpots[3*spotNum+1];
        double rad = loader->spotsRadius[spotNum];
        
        double sphereZ = sqrt(rad*rad - tx*tx - ty*ty);
        double R = sqrt(tx*tx + ty*ty);
        
        double phi = tx == 0.0 ? RNM : (atan2(ty, tx) * 180 / PI);
        
        double psi = R == 0.0 ? RNM : (atan2(sphereZ, R) * 180 / PI);
        
        double Vx, Vy, Vz;
        
        double locFac = 0.005;
        
        double cos_phi = R == 0.0 ? RNM : tx/R;
        double sin_phi = R == 0.0 ? RNM : ty/R;
        
        double cos_psi = R/rad;
        double sin_psi = sqrt(1-cos_psi*cos_psi);
        
        
        switch (loader->dim) {
            case 1:
                predVel[0] += V0 * locFac;
                corrVel[0] += V0 * locFac;
                
                break;
            case 2:
                psi = tx > 0.0 ? 90 : -90;
                
                Vx = V0 * cos(phi) * sin(psi);
                Vy = V0 * sin(phi);
                
                predVel[0] += Vx * locFac;
                corrVel[0] += Vx * locFac;
                predVel[1] += Vy * locFac;
                corrVel[1] += Vy * locFac;
                
                break;
            case 3:
                Vx = V0 * cos_phi * cos_psi;
                Vy = V0 * sin_phi * cos_psi;
                Vz = V0 * sin_psi;
                predVel[0] += Vx * locFac;
                corrVel[0] += Vx * locFac;
                predVel[1] += Vy * locFac;
                corrVel[1] += Vy * locFac;
                predVel[2] += Vz * locFac;
                corrVel[2] += Vz * locFac;
                
                break;
                
            default:
                throw runtime_error("dimension is not 1/2/3");
        }
        
        pusher->setParticleVelocity(m, 0, predVel[0]);
        pusher->setParticleVelocity(m, 3, corrVel[0]);
        
        pusher->setParticleVelocity(m, 1, predVel[1]);
        pusher->setParticleVelocity(m, 4, corrVel[1]);
        
        pusher->setParticleVelocity(m, 2, predVel[2]);
        pusher->setParticleVelocity(m, 5, corrVel[2]);
        
        //position of front from blast wave theory
        //
        // R = Eps0 (E/rho0)^0.2 t^0.4
        // dtR = 0.4 Eps0 (E/rho0)^0.2 t^-0.6
        // Eps0 ~ 1 geometrical factor
        // rho0 density of bg
        // E plasma energy = nT
        // t - time
        //
        // valid for R << (E/P0)^1/3
        // P0 - ambient gas pressure = 0.4*0.001 = 0.0004
        //
        // when
        // pressure driving shock wave
        // becomes comparable to ambient gas pressure P0
        // shock wave transforms in a sound wave
    }
    
    
    
    VectorVar** presNow = gridMgr->getVectorVariableOnG2(PRESSURE);
    VectorVar** densNow = gridMgr->getVectorVariableOnG2(DENSELEC);
    double G2shift = 0.5;
    // update electrons
    for (i = 1; i < xRes + 1; i++) {
        for (j = 1; j < yRes + 1; j++) {
            for (k = 1; k < zRes + 1; k++) {
                pos[0] = (i - G2shift) * dx + domainShiftX;
                pos[1] = (j - G2shift) * dy + domainShiftY;
                pos[2] = (k - G2shift) * dz + domainShiftZ;
                
                eleTemp(pos, T);
                
                if(T[0] <= 0.0){
                    continue;
                }
                idxOnG2 = IDX(i  , j  , k  , xResG2, yResG2, zResG2);
                
                currentDens = densNow[idxOnG2]->getValue()[0];
                double pres = currentDens*T[0];
                const double* pnow = presNow[idxOnG2]->getValue();
                double deltaP, rate, pres2set;
                rate = 0.0005;
                deltaP = pres;
                pres2set = deltaP*rate;
                gridMgr->addVectorVariableForNodeG2(idxOnG2, PRESSURE,     0, pres2set);
                gridMgr->addVectorVariableForNodeG2(idxOnG2, PRESSURE_AUX, 0, pres2set);
                gridMgr->addVectorVariableForNodeG2(idxOnG2, PRESSURE,     3, pres2set);
                gridMgr->addVectorVariableForNodeG2(idxOnG2, PRESSURE_AUX, 3, pres2set);
                gridMgr->addVectorVariableForNodeG2(idxOnG2, PRESSURE,     5, pres2set);
                gridMgr->addVectorVariableForNodeG2(idxOnG2, PRESSURE_AUX, 5, pres2set);
                
            }
        }
    }

    gridMgr->sendBoundary2Neighbor(PRESSURE);
    gridMgr->applyBC(PRESSURE);

    auto end_time = high_resolution_clock::now();
    string msg ="[LaserMockManager] accelerate() duration = "
                +to_string(duration_cast<milliseconds>(end_time - start_time).count())+" ms";
    logger->writeMsg(msg.c_str(), DEBUG);
    
}



