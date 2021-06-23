#include "ModelInitializer.hpp"



using namespace std;
using namespace chrono;




ModelInitializer::ModelInitializer(shared_ptr<Loader> load,
                                   shared_ptr<GridManager> grid,
                                   shared_ptr<Pusher> push):
                                   loader(move(load)), gridMng(move(grid)),
                                   pusher(move(push)){
    logger.reset(new Logger());
    
    initialize();
    logger->writeMsg("[ModelInitializer] create...OK", DEBUG);
}


void ModelInitializer::initialize()
{
    logger->writeMsg("[ModelInitializer] initialize...OK", DEBUG);
    
    totParticleNumberInDomain = 0;
    totalDensityInTHEbox = 0.0;
    
    switch (loader->runType) {
        case RESTART:
            readAllFromFile();
            break;
        case SCRATCH:
            initMagneticField();
            initVariablesonG2();// pressure density
            initParticles();// set coordinates velocities
            break;
            
        default:
            throw runtime_error("no runType");
    }
   
}


void ModelInitializer::readAllFromFile(){

    int rank ;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    
    int x, y, z;
    int i, j, k;
    int idx, idxG1;
    
    hid_t group, file, dataset, filespace;
    string fileName = loader->inputfile;
    
    file = H5Fopen(fileName.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    
    string groupname = "/vars";
    group = H5Gopen(file, groupname.c_str(), H5P_DEFAULT);
    
    int ptclNum = 0;
    
    int coresNum;
    MPI_Comm_size(MPI_COMM_WORLD, &coresNum);
    
    int* g2nodes  = new int[coresNum];
    int* g2offset = new int[coresNum];
   
    int* g1nodes  = new int[coresNum];
    int* g1offset = new int[coresNum];
    
    int* parts      = new int[coresNum];
    int* partOffset = new int[coresNum];

    hsize_t numOfCores = coresNum;
    hid_t memspaceAttr = H5Screate_simple(1,&numOfCores,NULL);

    readField(group, "g2_num"   , H5T_NATIVE_INT, file, memspaceAttr, g2nodes);
    readField(group, "g2_offset", H5T_NATIVE_INT, file, memspaceAttr, g2offset);

    readField(group, "g1_num"   , H5T_NATIVE_INT, file, memspaceAttr, g1nodes);
    readField(group, "g1_offset", H5T_NATIVE_INT, file, memspaceAttr, g1offset);
    
    readField(group, "parts_num"   , H5T_NATIVE_INT, file, memspaceAttr, parts);
    readField(group, "parts_offset", H5T_NATIVE_INT, file, memspaceAttr, partOffset);
    
    int numOfSpecies = loader->getNumberOfSpecies();
    pusher->initParticles(parts[rank], numOfSpecies);
    
    double* weights = new double[coresNum];
    
    for( int spn = 0; spn < numOfSpecies; spn++ ){
        pusher->setParticleMass4Type(spn, loader->getMass4species(spn));
        pusher->setParticleCharge4Type(spn, loader->getCharge4species(spn));
        pusher->setIfParticleTypeIsFrozen(spn, loader->getIfSpeciesFrozen(spn));
        
        readField(group, ("weight_"+to_string(spn)).c_str(), H5T_NATIVE_DOUBLE, file,
                  memspaceAttr, weights);
        pusher->setParticleWeight4Type(spn, weights[rank]);
    }
    delete[] weights;
    
    int totVarsOnG2 = gridMng->getVarsNumOnG2();
    hsize_t loc = g2nodes[rank];
    hsize_t memspace ;
    hsize_t locOffset = g2offset[rank];
    
    for ( int varN = 0; varN < totVarsOnG2; varN++ ){
        VectorVar** vars = gridMng->getVectorVariableOnG2(varN);
        int varSize = vars[0]->getSize();
        
        double* field = new double[g2nodes[rank]];
         for ( int dir = 0; dir < varSize; dir++ ){
             
            hid_t data = H5Dopen(group, ("g2_"+to_string(varN)+"_"+to_string(dir)).c_str(), H5P_DEFAULT);
            hid_t dataspace = H5Dget_space(data);
            H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, &locOffset, NULL, &loc, NULL);
            memspace  = H5Screate_simple(1, &loc , NULL);
            H5Dread(data, H5T_NATIVE_DOUBLE, memspace, dataspace, H5P_DEFAULT, field);

            for ( idx = 0; idx < g2nodes[rank]; idx++ ){
                 gridMng->setVectorVariableForNodeG2(idx, varN, dir, field[idx]);
            }
             
            H5Dclose(data);
            H5Sclose(memspace);
         }
        delete[] field;
    }
    

    locOffset = g1offset[rank];
    loc = g1nodes[rank];
    
    int totVarsOnG1 = gridMng->getVarsNumOnG1();
    
    for ( int varN=0; varN<totVarsOnG1; varN++){
        VectorVar** vars = gridMng->getVectorVariableOnG1(varN);
        int varSize = vars[0]->getSize();
        
        double* field = new double[g1nodes[rank]];
        for ( int dir=0; dir<varSize; dir++){
            
            hid_t data = H5Dopen(group, ("g1_"+to_string(varN)+"_"+to_string(dir)).c_str(), H5P_DEFAULT);
            hid_t dataspace = H5Dget_space(data);
            H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, &locOffset, NULL, &loc, NULL);
            memspace  = H5Screate_simple(1, &loc , NULL);
            H5Dread(data, H5T_NATIVE_DOUBLE, memspace, dataspace, H5P_DEFAULT, field);
  
            for ( idx=0; idx<g1nodes[rank]; idx++){
                gridMng->setVectorVariableForNodeG1(idx, varN, dir, field[idx]);
            }
            H5Dclose(data);
            H5Sclose(memspace);
        }
        delete[] field;
    }
    
    Particle** particles = pusher->getParticles();
    
    int totalPrtclNumber = pusher->getTotalParticleNumber();
    
    double* particlesFromFIle = new double[PARTICLES_SIZE*totalPrtclNumber];
    hid_t data = H5Dopen(group, "parts", H5P_DEFAULT);
    hid_t dataspace = H5Dget_space(data);
    
    loc = parts[rank]*PARTICLES_SIZE;
    locOffset = partOffset[rank]*PARTICLES_SIZE;
    
    H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, &locOffset, NULL, &loc, NULL);
    memspace  = H5Screate_simple(1, &loc , NULL);
    H5Dread(data, H5T_NATIVE_DOUBLE, memspace, dataspace, H5P_DEFAULT, particlesFromFIle);
   
    for ( idx=0; idx<totalPrtclNumber; idx++){
        particles[idx]->deserialize(particlesFromFIle, idx*PARTICLES_SIZE);
        
    }
    
    H5Dclose(data);
    H5Sclose(memspace);
    H5Sclose(memspaceAttr);
    H5Gclose(group);
    H5Fclose(file);
    
    delete[] particlesFromFIle;
    delete[] g2nodes;
    delete[] g2offset;
    delete[] g1nodes;
    delete[] g1offset;
    delete[] parts;
    delete[] partOffset;

}

void ModelInitializer::readField(hid_t group, string field,
                                 hid_t type, hid_t file,
                                 hid_t memspace, void * buf){
    
    hid_t data = H5Dopen(group, field.c_str(), H5P_DEFAULT);
    hid_t filespace = H5Dget_space(data);
    H5Dread(data, type, memspace, filespace, H5P_DEFAULT, buf);
    H5Dclose(data);
    
}


void ModelInitializer::initMagneticField(){
    logger->writeMsg("[ModelInitializer] init magnetic field...", DEBUG);
    int xRes = loader->resolution[0];
    int yRes = loader->resolution[1];
    int zRes = loader->resolution[2];
    
    double domainShiftX = loader->boxCoordinates[0][0];
    double domainShiftY = loader->boxCoordinates[1][0];
    double domainShiftZ = loader->boxCoordinates[2][0];
    
    double cellSizeX = loader->spatialSteps[0];
    double cellSizeY = loader->spatialSteps[1];
    double cellSizeZ = loader->spatialSteps[2];
    
    int xResG1 = xRes+1, yResG1 = yRes+1, zResG1 = zRes+1;
    double x,y,z;
    vector<double> bField;
    
    int i,j,k,idx;

    for ( i=0; i<xResG1; i++){
        for ( j=0; j<yResG1; j++){
            for ( k=0; k<zResG1; k++){
                idx = IDX(i,j,k,xResG1,yResG1,zResG1);
                
                x = i*cellSizeX + domainShiftX;
                y = j*cellSizeY + domainShiftY;
                z = k*cellSizeZ + domainShiftZ;
                
                bField = loader->getBfield(x, y, z);
                
                gridMng->setVectorVariableForNodeG1(idx, VectorVar(MAGNETIC    , bField));
                gridMng->setVectorVariableForNodeG1(idx, VectorVar(MAGNETIC_AUX, bField));

                
            }
        }
    }
    logger->writeMsg("[ModelInitializer] init magnetic field...OK", DEBUG);
    
}
                                         


void ModelInitializer::initVariablesonG2(){
    logger->writeMsg("[ModelInitializer] init G2...", DEBUG);
    int xRes = loader->resolution[0];
    int yRes = loader->resolution[1];
    int zRes = loader->resolution[2];
    
    double domainShiftX = loader->boxCoordinates[0][0];
    double domainShiftY = loader->boxCoordinates[1][0];
    double domainShiftZ = loader->boxCoordinates[2][0];
    
    double cellSizeX = loader->spatialSteps[0];
    double cellSizeY = loader->spatialSteps[1];
    double cellSizeZ = loader->spatialSteps[2];
    
    int xResG2 = xRes+2, yResG2 = yRes+2, zResG2 = zRes+2;
    double x,y,z;
    
    double G2shift = 0.5;
    int numOfSpecies = loader->getNumberOfSpecies();
    
    double densTot = 0.0, densTotLoc = 0.0, dens = 0.0, pres;
    int i,j,k,idx, spn;
    vector<double> densDomain;
    vector<double> locMIN4species;// need to provide ppc in the lowest density regions
    vector<double> locMAX4species;
    // need to provide ppc if lowest density is 0.0,
    // by default consider constant density everywhere
    
    prtcleWeight.reserve(numOfSpecies);
    
    for( spn = 0; spn < numOfSpecies; spn++ ){
        densDomain.push_back(0.0);
        locMIN4species.push_back(BIGN);//big number
        locMAX4species.push_back(0.0);
        prtcleWeight.push_back(0.0);
        npc.push_back(0);
    }
    
    int type2load = loader->prtclType2Load;
    
    // last element is reserved for loaded particles
    // -1 because there is no need to initialize zero number of loaded particles
    if ( loader->numOfSpots > 0 ) {
        numOfSpecies -= 1;
    }
    
    for( i = 0; i < xRes; i++ ){
        for( j = 0; j < yRes; j++ ){
            for( k = 0; k < zRes; k++ ){
                
                idx = IDX(i+1,j+1,k+1,xResG2,yResG2,zResG2);
                
                x = (i + G2shift)*cellSizeX + domainShiftX;
                y = (j + G2shift)*cellSizeY + domainShiftY;
                z = (k + G2shift)*cellSizeZ + domainShiftZ;
                
                densTotLoc = 0.0;
                for( spn = 0; spn < numOfSpecies; spn++ ){
                    dens = loader->getDensity(x, y, z, spn);
                    
                    // need only for particles initialization
                    // density is managed by HydroManager.cpp
                    // no need in boundary MPI communication
                    
                    if( loader->numOfSpots > 0 && spn == type2load){// replace particles by special type for loaded
                        double pres = loader->getElectronPressureProfile(x,y,z);
                        if( pres > 0.0 ){
                            gridMng->setVectorVariableForNodeG2(idx, gridMng->DENS_VEL(numOfSpecies), 0, dens);
                            densDomain[numOfSpecies] += dens;
                        }else{
                            gridMng->setVectorVariableForNodeG2(idx, gridMng->DENS_VEL(spn), 0, dens);
                            densDomain[spn] += dens;
                        }
                    }else{
                        gridMng->setVectorVariableForNodeG2(idx, gridMng->DENS_VEL(spn), 0, dens);
                        densDomain[spn] += dens;
                    }
                    
                    densTotLoc += dens;
                    if( dens < locMIN4species[spn] ){
                        locMIN4species[spn] = dens;
                    }
                    if( dens > locMAX4species[spn] ){
                        locMAX4species[spn] = dens;
                    }
                }
                densTot += densTotLoc;
                
            }
        }
    }
    
    
    int rank ;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    
    double localMinimumDens, localMaximumDens;
    int nonemptyCellLoc, nonemptyCellGlob;
    int partclNumPerDomain;
    double globalMaximumDens;
    
    if ( loader->numOfSpots > 0 ) {
        numOfSpecies += 1;// set back normal value
    }
    
    
    int sum = 0;
    for( spn = 0; spn < numOfSpecies; spn++ ){
        
        // need this routine for non uniform plasma distribution
        localMinimumDens = locMIN4species[spn];
        MPI_Allreduce(&localMinimumDens, &globalMinimumDens, 1, MPI_DOUBLE_INT, MPI_MINLOC, MPI_COMM_WORLD);
        localMaximumDens = locMAX4species[spn];
        MPI_Allreduce(&localMaximumDens, &globalMaximumDens, 1, MPI_DOUBLE_INT, MPI_MAXLOC, MPI_COMM_WORLD);
        
        if( globalMinimumDens <= 0.0 ){
            globalMinimumDens = localMaximumDens > 0.0 ? localMaximumDens : globalMaximumDens;
        }
        
        if( globalMinimumDens <= loader->minimumDens2ResolvePPC ){
            globalMinimumDens = loader->minimumDens2ResolvePPC;
        }
        
        int ppc = loader->getPPC4species(spn);
        prtcleWeight[spn] = globalMinimumDens/ppc;
        
        if ( (loader->numOfSpots > 0) && (spn == (numOfSpecies-1)) ) {
            
            int type2load = loader->prtclType2Load;
            int ppc4loadedType = loader->getPPC4species(type2load);
            
            const double weightDecreaseFactor = double(ppc4loadedType)/double(ppc);
            
            double loadWeight = weightDecreaseFactor*prtcleWeight[type2load];
            
            prtcleWeight[spn] = loadWeight;
        }
        
        partclNumPerDomain = int(densDomain[spn]/prtcleWeight[spn]);
        npc[spn] = partclNumPerDomain;
        sum += partclNumPerDomain;
        
        // additional check for extraordinary situations
        double localWeight = prtcleWeight[spn];
        double globWeightMax;
        double globWeightMin;
        
        MPI_Allreduce(&localWeight, &globWeightMin, 1, MPI_DOUBLE_INT, MPI_MINLOC, MPI_COMM_WORLD);
        MPI_Allreduce(&localWeight, &globWeightMax, 1, MPI_DOUBLE_INT, MPI_MAXLOC, MPI_COMM_WORLD);
        
        if( globWeightMax != globWeightMin ){
            throw runtime_error(" Different particle weight on different domains!");
        }

        
        logger->writeMsg(("[ModelInitializer] particle weight for species "+to_string(spn+1)
                          +" = "+to_string(prtcleWeight[spn])
                          +"\n"+string( 17, ' ' )+" required particles number for domain "+to_string(rank)
                          +" =  "+to_string(partclNumPerDomain)
                          +"\n"+string( 17, ' ' )+" globalMinimumDens  = "+to_string(globalMinimumDens)
                          +"\n"+string( 17, ' ' )+" globalMaximumDens  = "+to_string(globalMaximumDens)
                          +"\n"+string( 17, ' ' )+" localMinimumDens  = "+to_string(localMinimumDens)
                          +"\n"+string( 17, ' ' )+" localMaximumDens  = "+to_string(localMaximumDens)
                          +"\n"+string( 17, ' ' )+" localMaximumDens  = "+to_string(localMaximumDens)
                          +"\n"+string( 17, ' ' )+" densDomain  = "+to_string(densDomain[spn])
                          +"\n"+string( 17, ' ' )+" ppc  = "+to_string(ppc)).c_str(),  DEBUG);
    }
    totParticleNumberInDomain = sum;
    MPI_Allreduce(&densTot, &totalDensityInTHEbox, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    string msg00 ="[ModelInitializer] total density  = "+to_string(totalDensityInTHEbox)+" totParticleNumberInDomain = "+to_string(totParticleNumberInDomain);
    logger->writeMsg(msg00.c_str(),  DEBUG);
    
    logger->writeMsg("[ModelInitializer] init G2...OK", DEBUG);
    
}



void ModelInitializer::initParticles(){
    logger->writeMsg("[ModelInitializer] init particles...", DEBUG);
    double domainShiftX = loader->boxCoordinates[0][0];
    double domainShiftY = loader->boxCoordinates[1][0];
    double domainShiftZ = loader->boxCoordinates[2][0];
    
    int xRes = loader->resolution[0];
    int yRes = loader->resolution[1];
    int zRes = loader->resolution[2];
    
    double Lx = loader->boxSizes[0], Ly = loader->boxSizes[1], Lz = loader->boxSizes[2];
    
    double cellSizeX = loader->spatialSteps[0];
    double cellSizeY = loader->spatialSteps[1];
    double cellSizeZ = loader->spatialSteps[2];
    
    int numOfSpecies = loader->getNumberOfSpecies();
    
    int particle_idx=0, idxOnG2, idx, ptclIDX, spn;
    double pos[3], vpb[3];
    vector<double> vel;
    double r1, r2;
    
    
    if( totParticleNumberInDomain <= 0 ){
        throw runtime_error("Check initialization, particles number in domain = "
                            +to_string(totParticleNumberInDomain));
    }
    
    //first create particles with default parameters
    pusher->initParticles(totParticleNumberInDomain, numOfSpecies);
    
    double totalPrtclNumber = pusher->getTotalParticleNumber();
    
    string msg001a1 ="[ModelInitializer] have to init totParticleNumberInDomain  = "
                        +to_string(totParticleNumberInDomain)
                        +"\n"+string( 19, ' ' )+"  globalMinimumDens = "+to_string(globalMinimumDens)
                        +"\n"+string( 19, ' ' )+"  totalPrtclNumber from pusher = "+to_string(totalPrtclNumber);
    
    logger->writeMsg(msg001a1.c_str(),  DEBUG);
    
    map<int, VectorVar**> dens;
    
    for( spn = 0; spn < numOfSpecies; spn++ ){
        
        if ( (loader->numOfSpots > 0) && (spn == (numOfSpecies-1)) ) {
            int type2load = loader->prtclType2Load;
            double loadMass = pusher->getParticleMass4Type(type2load);
            double loadCharge = pusher->getParticleCharge4Type(type2load);
            // last element is reserved for loaded particles
            pusher->setParticleMass4Type(spn, loadMass);
            pusher->setParticleCharge4Type(spn, loadCharge);
        }else{
            pusher->setParticleMass4Type(spn, loader->getMass4species(spn));
            pusher->setParticleCharge4Type(spn, loader->getCharge4species(spn));
            pusher->setIfParticleTypeIsFrozen(spn, loader->getIfSpeciesFrozen(spn));
        }
        
        pusher->setParticleWeight4Type(spn, prtcleWeight[spn]);
        dens[spn] = gridMng->getVectorVariableOnG2(gridMng->DENS_VEL(spn));
    }
    
   
    vector<int> ppcValues;
    for( spn = 0; spn < numOfSpecies; spn++ ){
        int ppc = loader->getPPC4species(spn);
        ppcValues.push_back(ppc);
    }
    
    double requiredPrtclNum;
    
    int i, j, k;
    for( i = 0; i < xRes; i++){
        for( j = 0; j < yRes; j++) {
            for( k = 0; k < zRes; k++){
                
                idxOnG2 = IDX(i+1, j+1, k+1, xRes+2, yRes+2, zRes+2);
                
                for( spn = 0; spn < numOfSpecies; spn++ ){
                    
                    requiredPrtclNum = int(dens[spn][idxOnG2]->getValue()[0]/
                                           pusher->getParticleWeight4Type(spn));
                    
                    
                    for (ptclIDX=0; ptclIDX < requiredPrtclNum; ptclIDX++){
                        
                        double ran1, ran2, ran3;
                        
                        if( ppcValues[spn] == 1 ){
                            ran1 = 0.5;// put at the cell center
                            ran2 = 0.5;
                            ran3 = 0.5;
                        }else{
                            ran1 = RNM;
                            ran2 = RNM;
                            ran3 = RNM;
                        }

                        ran1 = ran1 == 1.0 ? 1-EPS4 : ran1;
                        ran2 = ran2 == 1.0 ? 1-EPS4 : ran2;
                        ran3 = ran3 == 1.0 ? 1-EPS4 : ran3;
                        
                        pos[0] = (i + ran1) * cellSizeX + domainShiftX;
                        pos[1] = (j + ran2) * cellSizeY + domainShiftY;
                        pos[2] = (k + ran3) * cellSizeZ + domainShiftZ;
                        
                        if( particle_idx >= totalPrtclNumber ){
                            string msg0011 ="[ModelInitializer] idx is out of preinitialized particles number: particle_idx  = "+to_string(particle_idx);
                            logger->writeMsg(msg0011.c_str(),  DEBUG);
                        }
                        
                        if( pos[0] < 0.0 || pos[1] < 0.0 || pos[2] < 0.0 ){
                            throw runtime_error("position pos[0] = "+to_string(pos[0]));
                        }
                        double pos2Save[6] = {pos[0], pos[1], pos[2], pos[0], pos[1], pos[2]};
                        pusher->setParticlePosition(particle_idx, pos2Save);
                        
                        pusher->setParticleType(particle_idx,  spn);

                        vel = loader->getVelocity(pos[0], pos[1], pos[2], spn);
                    
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
                        
                        double vel2Save[6] = {vpb[0], vpb[1], vpb[2], vpb[0], vpb[1], vpb[2]};
                        pusher->setParticleVelocity(particle_idx, vel2Save);
                        
                        particle_idx++;
                    }
                }
            }
        }
    }
    pusher->setTotalParticleNumber(particle_idx);

    logger->writeMsg("[ModelInitializer] init particles...OK", DEBUG);
}


ModelInitializer::~ModelInitializer(){
    finilize();
}

void ModelInitializer::finilize(){
    
}
