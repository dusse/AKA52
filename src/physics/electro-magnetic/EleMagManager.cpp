#include "EleMagManager.hpp"

using namespace std;
using namespace chrono;


EleMagManager::EleMagManager(std::shared_ptr<Loader> ldr,
                             std::shared_ptr<GridManager> gridMnr):loader(move(ldr)), gridMgr(move(gridMnr)){
    logger.reset(new Logger());
    initialize();
    logger->writeMsg("[EleMagManager] create...OK", DEBUG);
}

void EleMagManager::initialize(){
    
    int xRes = loader->resolution[0],
    yRes = loader->resolution[1],
    zRes = loader->resolution[2];
    int xResG2 = xRes+2, yResG2 = yRes+2, zResG2 = zRes+2;
    
    int nG2= xResG2*yResG2*zResG2;
    
    //magnetic field is defined by ModelInitializer
    
    calculateBhalf(PREDICTOR);
    calculateJhalf(PREDICTOR);
    
    VectorVar** current = gridMgr->getVectorVariableOnG2(CURRENT);
    int ijkG2, h;
    // CURRENT_AUX is smoothed CURRENT, normally CURRENT_AUX is set in ClosureManager.cpp
    for(ijkG2=0; ijkG2<nG2; ijkG2++){
        for (h = 0; h < 3; h++) {
            gridMgr->setVectorVariableForNodeG2(ijkG2, CURRENT_AUX, h,
                                                current[ijkG2]->getValue()[h]);
        }
    }
    
    calculateEnext(PREDICTOR);
    calculateBnext();
    calculateJnext();
    
    calculateBhalf(CORRECTOR);
    calculateJhalf(CORRECTOR);
    
    for(ijkG2=0; ijkG2<nG2; ijkG2++){
        for (h = 0; h < 3; h++) {
            gridMgr->setVectorVariableForNodeG2(ijkG2, CURRENT_AUX, h,
                                                current[ijkG2]->getValue()[h]);
        }
    }
    
    calculateEnext(CORRECTOR);
    calculateBnext();
    calculateJnext();

}


void EleMagManager::calculateMagneticField(int magField2use, int eleField2use, int magField2save){
    
    double dtx, dty, dtz;
    double ts = loader->getTimeStep();
    
    double dx = loader->spatialSteps[0];
    double dy = loader->spatialSteps[1];
    double dz = loader->spatialSteps[2];
    
    int xSize = loader->resolution[0];
    int ySize = loader->resolution[1];
    int zSize = loader->resolution[2];
    
    dtx = 0.125*ts/dx;
    dty = 0.125*ts/dy;
    dtz = 0.125*ts/dz;
    
    int totalG1Nodes = (xSize+1)*(ySize+1)*(zSize+1);

    const int * ngborsXder = gridMgr->getNeibors4G2spatialDerX();
    const int * ngborsYder = gridMgr->getNeibors4G2spatialDerY();
    const int * ngborsZder = gridMgr->getNeibors4G2spatialDerZ();
    
    VectorVar** eField     = gridMgr->getVectorVariableOnG2(eleField2use);
    VectorVar** bFieldPrev = gridMgr->getVectorVariableOnG1(magField2use);

    double Exdy, Exdz, Eydx, Eydz, Ezdx, Ezdy;
    int left, rigt;
    double bFieldX, bFieldY, bFieldZ;
    const double* B0;

    for(int idx = 0; idx<totalG1Nodes; idx++){
        
        Exdy = 0; Exdz= 0; Eydx= 0; Eydz= 0; Ezdx= 0; Ezdy = 0;
        for (int pairNum = 0; pairNum<4; pairNum++){
            left = ngborsXder[8*idx+2*pairNum+0];
            rigt = ngborsXder[8*idx+2*pairNum+1];// index ijk is saved with 1
            Eydx += eField[left]->getValue()[1] - eField[rigt]->getValue()[1];
            Ezdx += eField[left]->getValue()[2] - eField[rigt]->getValue()[2];
                
            left = ngborsYder[8*idx+2*pairNum+0];
            rigt = ngborsYder[8*idx+2*pairNum+1];
            Exdy += eField[left]->getValue()[0] - eField[rigt]->getValue()[0];
            Ezdy += eField[left]->getValue()[2] - eField[rigt]->getValue()[2];
            
            left = ngborsZder[8*idx+2*pairNum+0];
            rigt = ngborsZder[8*idx+2*pairNum+1];
            Exdz += eField[left]->getValue()[0] - eField[rigt]->getValue()[0];
            Eydz += eField[left]->getValue()[1] - eField[rigt]->getValue()[1];
        }
        B0 = bFieldPrev[idx]->getValue();
        bFieldX = B0[0] + Eydz*dtz - Ezdy*dty;
        bFieldY = B0[1] + Ezdx*dtx - Exdz*dtz;
        bFieldZ = B0[2] + Exdy*dty - Eydx*dtx;
        
        gridMgr->setVectorVariableForNodeG1(idx, magField2save, 0, bFieldX);
        gridMgr->setVectorVariableForNodeG1(idx, magField2save, 1, bFieldY);
        gridMgr->setVectorVariableForNodeG1(idx, magField2save, 2, bFieldZ);

    }
    
}


void EleMagManager::calculateCurrent(int magField2use, int current2save){
    
    int idxG1, idxG2;
    double dx = loader->spatialSteps[0];
    double dy = loader->spatialSteps[1];
    double dz = loader->spatialSteps[2];
    
    int xSize = loader->resolution[0];
    int ySize = loader->resolution[1];
    int zSize = loader->resolution[2];
    
    dx = 1.0/dx;
    dy = 1.0/dy;
    dz = 1.0/dz;
    
    const int* ngborsXder = gridMgr->getNeibors4G1spatialDerX();
    const int* ngborsYder = gridMgr->getNeibors4G1spatialDerY();
    const int* ngborsZder = gridMgr->getNeibors4G1spatialDerZ();
    
    VectorVar** bField = gridMgr->getVectorVariableOnG1(magField2use);
    
    double Bxdy, Bxdz, Bydx, Bydz, Bzdx, Bzdy;
    int left, rigt;
    double currentX, currentY, currentZ;
    
    int i,j,k;
    for ( i=1; i<xSize+1; i++){
        for ( j=1; j<ySize+1; j++){
            for ( k=1; k<zSize+1; k++){
                
                idxG1 = IDX(i,j,k,xSize+1,ySize+1,zSize+1);
                
                Bxdy = 0; Bxdz = 0;
                Bydx = 0; Bydz = 0;
                Bzdx = 0; Bzdy = 0;
        
                for (int pairNum = 0; pairNum<4; pairNum++){
            
                    left = ngborsXder[8*idxG1+2*pairNum+0];
                    rigt = ngborsXder[8*idxG1+2*pairNum+1];// index ijk is saved in +1
                    
                    Bydx += bField[rigt]->getValue()[1]
                          - bField[left]->getValue()[1];
                    
                    Bzdx += bField[rigt]->getValue()[2]
                          - bField[left]->getValue()[2];
            
                    left = ngborsYder[8*idxG1+2*pairNum+0];
                    rigt = ngborsYder[8*idxG1+2*pairNum+1];
            
                    Bxdy += bField[rigt]->getValue()[0]
                          - bField[left]->getValue()[0];
                    
                    Bzdy += bField[rigt]->getValue()[2]
                          - bField[left]->getValue()[2];
            
                    left = ngborsZder[8*idxG1+2*pairNum+0];
                    rigt = ngborsZder[8*idxG1+2*pairNum+1];
            
                    Bxdz += bField[rigt]->getValue()[0]
                          - bField[left]->getValue()[0];
                    
                    Bydz += bField[rigt]->getValue()[1]
                          - bField[left]->getValue()[1];
            
                }
                currentX = 0.25*(Bzdy*dy - Bydz*dz);
                currentY = 0.25*(Bxdz*dz - Bzdx*dx);
                currentZ = 0.25*(Bydx*dx - Bxdy*dy);
        
                idxG2 = IDX(i,j,k,xSize+2,ySize+2,zSize+2);
        
                gridMgr->setVectorVariableForNodeG2(idxG2, current2save, 0, currentX);
                gridMgr->setVectorVariableForNodeG2(idxG2, current2save, 1, currentY);
                gridMgr->setVectorVariableForNodeG2(idxG2, current2save, 2, currentZ);
            }
        }
    }
    gridMgr->sendBoundary2Neighbor(current2save);
    
    gridMgr->applyBC(current2save);

}

void EleMagManager::calculateBhalf(int phase){
    auto start_time = high_resolution_clock::now();
    
    string msg;
    
    int magField2use;
    int magFeld2save;
    int eleField2use = ELECTRIC;
    
    switch (phase)
    {
        case PREDICTOR:
            magField2use = MAGNETIC;
            magFeld2save = MAGNETIC_AUX;
            break;
            
        case CORRECTOR:
            magField2use = MAGNETIC;
            magFeld2save = MAGNETIC;
            break;
            
        default :
            msg ="[EleMagManager] calculateBhalf() no such a phase = "+to_string(phase);
            logger->writeMsg(msg.c_str(), CRITICAL);
            throw runtime_error("no phase");
    }
    
            
    calculateMagneticField(magField2use, eleField2use, magFeld2save);
    
    auto end_time = high_resolution_clock::now();
    msg ="[EleMagManager] calculateBhalf() duration = "
            +to_string(duration_cast<milliseconds>(end_time - start_time).count())+" ms";
    logger->writeMsg(msg.c_str(), DEBUG);
}


void EleMagManager::calculateBnext(){
    auto start_time = high_resolution_clock::now();
    
    int magField2use = MAGNETIC_AUX;
    int magFeld2save = MAGNETIC;
    int eleField2use = ELECTRIC;
    
    calculateMagneticField(magField2use, eleField2use, magFeld2save);
    
    auto end_time = high_resolution_clock::now();
    string msg ="[EleMagManager] calculateBnext() duration = "
                +to_string(duration_cast<milliseconds>(end_time - start_time).count())+" ms";
    logger->writeMsg(msg.c_str(), DEBUG);
}


void EleMagManager::calculateJhalf(int phase){
    auto start_time = high_resolution_clock::now();
    
    int magField2use = MAGNETIC;
    int current2save = CURRENT;
    
    calculateCurrent(magField2use, current2save);
    
    auto end_time = high_resolution_clock::now();
    string msg ="[EleMagManager] calculateJhalf() duration = "
                +to_string(duration_cast<milliseconds>(end_time - start_time).count())+" ms";
    logger->writeMsg(msg.c_str(), DEBUG);
}

void EleMagManager::calculateJnext(){
    auto start_time = high_resolution_clock::now();
    
    int magField2use = MAGNETIC;
    int current2save = CURRENT;
    
    calculateCurrent(magField2use, current2save);
    
    auto end_time = high_resolution_clock::now();
    string msg ="[EleMagManager] calculateJnext() duration = "
                +to_string(duration_cast<milliseconds>(end_time - start_time).count())+" ms";
    logger->writeMsg(msg.c_str(), DEBUG);
}



void EleMagManager::calculateEnext(int phase){
    auto start_time = high_resolution_clock::now();
    
    double dx = loader->spatialSteps[0];
    double dy = loader->spatialSteps[1];
    double dz = loader->spatialSteps[2];
    
    int xSize = loader->resolution[0];
    int ySize = loader->resolution[1];
    int zSize = loader->resolution[2];
    
    int totG2 = (xSize+2)*(ySize+2)*(zSize+2);

    int magField2use, eleField2save;
    string msg;

    switch (phase)
    {
        case PREDICTOR:
            magField2use  = MAGNETIC_AUX;
            eleField2save = ELECTRIC_AUX;
            break;
            
        case CORRECTOR:
            magField2use  = MAGNETIC;
            eleField2save = ELECTRIC;
            break;
            
        default :
            msg ="[EleMagManager] calculateEnext() no such a phase = "+to_string(phase);
            logger->writeMsg(msg.c_str(), CRITICAL);
            throw runtime_error("no phase");
    }
    
    VectorVar** bField   = gridMgr->getVectorVariableOnG1(magField2use);
    VectorVar** presEle  = gridMgr->getVectorVariableOnG2(PRESSURE);
    VectorVar** current  = gridMgr->getVectorVariableOnG2(CURRENT_AUX);
    VectorVar** density  = gridMgr->getVectorVariableOnG2(DENSELEC);
    VectorVar** velocity = gridMgr->getVectorVariableOnG2(VELOCION);
    
    int compIDX[3][3]={{0, 1, 2},{1,3,4},{2,4,5}};
    double WEIHTS[9] = {0.125, 0.125, 0.125, 0.125, 0.250, 0.250, 0.250, 0.250, 0.500};
    double deltas[3] = {0.25/dx, 0.25/dy, 0.25/dz};
    double deltasLap[3] = {1.0/dx/dx, 1.0/dy/dy, 1.0/dz/dz};
    double wei;
    int left, rigt;
    
    const int* negbors4PresX = gridMgr->getNghbd4DivPeInXonG2();
    const int* negbors4PresY = gridMgr->getNghbd4DivPeInYonG2();
    const int* negbors4PresZ = gridMgr->getNghbd4DivPeInZonG2();
    
    const int* negbors4LaplacX = gridMgr->getNghbd4LaplacInXonG2();
    const int* negbors4LaplacY = gridMgr->getNghbd4LaplacInYonG2();
    const int* negbors4LaplacZ = gridMgr->getNghbd4LaplacInZonG2();
    
    const int* neighbourhood = gridMgr->getNeighbourhoodOnG1();
    
    int idxG1, idxG2, idxNeigbor, coord, neighbour, curPcomp;
    
    double locB[3], locE[3], lapJ[3], divP[3], dens;
    const double hypVis = loader->hyperviscosity;
    double dL, dP;
    
    
    int i,j,k;
    for ( i=1; i<xSize+1; i++){
        for ( j=1; j<ySize+1; j++){
            for ( k=1; k<zSize+1; k++){
                
                idxG1 = IDX(i,j,k,xSize+1,ySize+1,zSize+1);
                idxG2 = IDX(i,j,k,xSize+2,ySize+2,zSize+2);
                
        
                divP[0] = 0.0, divP[1] = 0.0, divP[2] = 0.0;
                lapJ[0] = 0.0; lapJ[1] = 0.0; lapJ[2] = 0.0;
                locB[0] = 0.0; locB[1] = 0.0; locB[2] = 0.0;
                locE[0] = 0.0; locE[1] = 0.0; locE[2] = 0.0;
        
                for (coord=0; coord<3; coord++){
            
                    for (neighbour=0; neighbour<8; neighbour++){
                        idxNeigbor = neighbourhood[8*idxG1+neighbour];
                        locB[coord] += 0.125*bField[idxNeigbor]->getValue()[coord];
                    }
            
                    dP = deltas[coord];
                    for (neighbour=0; neighbour<9; neighbour++){
                
                        wei  = WEIHTS[neighbour];
                    
                        curPcomp = compIDX[0][coord];
                        left = negbors4PresX[18*idxG2+2*neighbour+0];
                        rigt = negbors4PresX[18*idxG2+2*neighbour+1];
                        divP[0] += (presEle[left]->getValue()[curPcomp]
                                    -presEle[rigt]->getValue()[curPcomp])*wei*dP;
                        
                        curPcomp = compIDX[1][coord];
                        left = negbors4PresY[18*idxG2+2*neighbour+0];
                        rigt = negbors4PresY[18*idxG2+2*neighbour+1];
                        divP[1] += (presEle[left]->getValue()[curPcomp]
                                    -presEle[rigt]->getValue()[curPcomp])*wei*dP;
                
                        curPcomp = compIDX[2][coord];
                        left = negbors4PresZ[18*idxG2+2*neighbour+0];
                        rigt = negbors4PresZ[18*idxG2+2*neighbour+1];
                        divP[2] += (presEle[left]->getValue()[curPcomp]
                                    -presEle[rigt]->getValue()[curPcomp])*wei*dP;
                    }
            
                    dL = deltasLap[coord];
                    for (neighbour=0; neighbour<2; neighbour++){
                        left = negbors4LaplacX[4*idxG2+2*neighbour+0];
                        rigt = negbors4LaplacX[4*idxG2+2*neighbour+1];
                        lapJ[0] += (current[left]->getValue()[coord]
                                    -current[rigt]->getValue()[coord])*dL;
                        
                        left = negbors4LaplacY[4*idxG2+2*neighbour+0];
                        rigt = negbors4LaplacY[4*idxG2+2*neighbour+1];
                        lapJ[1] += (current[left]->getValue()[coord]
                                    -current[rigt]->getValue()[coord])*dL;
                
                        left = negbors4LaplacZ[4*idxG2+2*neighbour+0];
                        rigt = negbors4LaplacZ[4*idxG2+2*neighbour+1];
                        lapJ[2] += (current[left]->getValue()[coord]
                                    -current[rigt]->getValue()[coord])*dL;
                    }
                }

                dens = density[idxG2]->getValue()[0];
                const double* velI = velocity[idxG2]->getValue();
                const double* J    = current[idxG2]->getValue();
        
                locE[0] = - (velI[1]*locB[2] - velI[2]*locB[1])
                          + (   J[1]*locB[2] -    J[2]*locB[1])/dens
                          - divP[0]/dens
                          - hypVis*lapJ[0];
                
                locE[1] = - (velI[2]*locB[0] - velI[0]*locB[2])
                          + (   J[2]*locB[0] -    J[0]*locB[2])/dens
                          - divP[1]/dens
                          - hypVis*lapJ[1];
                
                locE[2] = - (velI[0]*locB[1] - velI[1]*locB[0])
                          + (   J[0]*locB[1] -    J[1]*locB[0])/dens
                          - divP[2]/dens
                          - hypVis*lapJ[2];

#ifdef LOG
                for(int coord=0; coord < 3; coord++){
                    if( std::isnan(locE[coord])){
                        write2Log(idxG2, i, j, k, velI, locB, divP, lapJ, J, dens);
                    }
                }
#endif
                
                gridMgr->setVectorVariableForNodeG2(idxG2, eleField2save, 0, locE[0]);
                gridMgr->setVectorVariableForNodeG2(idxG2, eleField2save, 1, locE[1]);
                gridMgr->setVectorVariableForNodeG2(idxG2, eleField2save, 2, locE[2]);
                
                
            }
        }
    }
 
    gridMgr->sendBoundary2Neighbor(eleField2save);
    
    gridMgr->applyBC(eleField2save);
    
    VectorVar** eField    = gridMgr->getVectorVariableOnG2(ELECTRIC);
    VectorVar** eFieldAux = gridMgr->getVectorVariableOnG2(ELECTRIC_AUX);
    
    
    locE[0] = 0.0; locE[1] = 0.0; locE[2] = 0.0;
    switch (phase){
            case PREDICTOR:
                for(int idx = 0; idx<totG2; idx++){
                    for (coord=0; coord<3; coord++){
                        locE[coord] =
                        - eField[idx]->getValue()[coord]+2.0*eFieldAux[idx]->getValue()[coord];
                        
                        gridMgr->setVectorVariableForNodeG2(idx, ELECTRIC, coord, locE[coord]);
                    }
                }
            break;
            case CORRECTOR:
                for(int idx = 0; idx<totG2; idx++){
                    for (coord=0; coord<3; coord++){
                        locE[coord] =
                        0.5*(eField[idx]->getValue()[coord]+eFieldAux[idx]->getValue()[coord]);
                        
                        gridMgr->setVectorVariableForNodeG2(idx, ELECTRIC, coord, locE[coord]);
                    }
                }
            break;
            default :
                throw runtime_error("no phase");
    }
    
    
    auto end_time = high_resolution_clock::now();
    msg ="[EleMagManager] calculateEnext() duration = "
            +to_string(duration_cast<milliseconds>(end_time - start_time).count())+" ms";
    logger->writeMsg(msg.c_str(), DEBUG);
}



void EleMagManager::write2Log(int idxG2, int i, int j, int k, const double* velI,
                              double* locB, double* divP, double* lapJ, const double* J, double dens){
    logger->writeMsg(("[EleMagManager] idxG2 =  "+to_string(idxG2)
                      +"\n     i  = "+to_string(i)
                      +"\n     j  = "+to_string(j)
                      +"\n     k  = "+to_string(k)
                      +"\n     Vion_x  = "+to_string(velI[0])
                      +"\n     Vion_y  = "+to_string(velI[1])
                      +"\n     Vion_z  = "+to_string(velI[2])
                      +"\n     Bx  = "+to_string(locB[0])
                      +"\n     By  = "+to_string(locB[1])
                      +"\n     Bz  = "+to_string(locB[2])
                      +"\n     divPx = "+to_string(divP[0])
                      +"\n     divPy = "+to_string(divP[1])
                      +"\n     divPz = "+to_string(divP[2])
                      +"\n     lapJx = "+to_string(lapJ[0])
                      +"\n     lapJy = "+to_string(lapJ[1])
                      +"\n     lapJz = "+to_string(lapJ[2])
                      +"\n     Jx = "+to_string(J[0])
                      +"\n     Jx = "+to_string(J[1])
                      +"\n     Jz = "+to_string(J[2])
                      +"\n     Ne  = "+to_string(dens)
                      ).c_str(), CRITICAL);
}





