#include "BoundaryManager.hpp"

using namespace std;
using namespace chrono;

BoundaryManager::BoundaryManager(shared_ptr<Loader> ldr):loader(move(ldr)){
    logger.reset(new Logger());
    initialize();
    string msg ="[BoundaryManager] init...OK";
    logger->writeMsg(msg.c_str(), DEBUG);
}



void BoundaryManager::initialize(){
    logger->writeMsg("[BoundaryManager] initialize() ...", DEBUG);
    leavingParticles.reserve(NUM_OF_LEAVING_PACTICLES);
    for( int t = 0; t < 27; t++ ){
        domain2send[t] = 0;
    }
    logger->writeMsg("[BoundaryManager] initialize() ...OK", DEBUG);
}



int BoundaryManager::isPtclOutOfDomain(double pos[3]){
    
    double domainShiftX = loader->boxCoordinates[0][0];
    double domainShiftY = loader->boxCoordinates[1][0];
    double domainShiftZ = loader->boxCoordinates[2][0];
    double domainXSize = loader->resolution[0];
    double domainYSize = loader->resolution[1];
    double domainZSize = loader->resolution[2];
    double dx = loader->spatialSteps[0];
    double dy = loader->spatialSteps[1];
    double dz = loader->spatialSteps[2];
    
    double x, y, z;
    int a, b, c, t;
    
    x = (pos[0] - domainShiftX)/domainXSize/dx;
    y = (pos[1] - domainShiftY)/domainYSize/dy;
    z = (pos[2] - domainShiftZ)/domainZSize/dz;
    
    switch (loader->dim) {
        case 1:
            a = (int)floor(x);
            b = 0;
            c = 0;
            break;
        case 2:
            a = (int)floor(x);
            b = (int)floor(y);
            c = 0;
            break;
        case 3:
            a = (int)floor(x);
            b = (int)floor(y);
            c = (int)floor(z);
            break;
            
        default:
            throw runtime_error("dimension is not 1/2/3");
    }
    
    t = (1+c)+3*((1+b)+3*(1+a));
    

    if ( t < 0 || t > 26 || std::isnan(x) || std::isnan(y) || std::isnan(z) ){
        
        t = IN; // keep particle on domain, pusher will remove it
        
        #ifdef HEAVYLOG
        string pb ="[BoundaryManager] problem particle destination:  t = "+to_string(t)
        +"\n      x = "+to_string(x)
        +"\n      y = "+to_string(y)
        +"\n      z = "+to_string(z)
        +"\n      domainShiftX = "+to_string(domainShiftX)
        +"\n      domainShiftY = "+to_string(domainShiftY)
        +"\n      domainShiftZ = "+to_string(domainShiftZ)
        +"\n      a = "+to_string(a)
        +"\n      b = "+to_string(b)
        +"\n      c = "+to_string(c);
        logger->writeMsg(pb.c_str(), CRITICAL);
        #endif
    }

    //left = 4 right = 22 top = 16 bottom = 10 back = 12 front = 14
    
    return (t == 13) ? IN : t;
}


void BoundaryManager::applyBC(Particle** particles,
                              vector<shared_ptr<Particle>> &particles2add,
                              int phase){
    auto start_time = high_resolution_clock::now();
    
    logger->writeMsg(("[BoundaryManager] applyBC .. leavingParticles = "
                     +to_string(leavingParticles.size())).c_str(), DEBUG);
    
    double domainShiftX = loader->boxCoordinates[0][0];
    double domainShiftY = loader->boxCoordinates[1][0];
    double domainShiftZ = loader->boxCoordinates[2][0];
    double domainXSize = loader->resolution[0];
    double domainYSize = loader->resolution[1];
    double domainZSize = loader->resolution[2];
    double dx = loader->spatialSteps[0];
    double dy = loader->spatialSteps[1];
    double dz = loader->spatialSteps[2];
    
    double Lx = loader->boxSizes[0],
           Ly = loader->boxSizes[1],
           Lz = loader->boxSizes[2];
    
    int idx, t, a, b, c;
    
    double xloc, yloc, zloc;
    double x0, y0, z0;
    
    double *sendBuf[27];
    double *recvBuf[27];
    int partcls2send[27];
    int partcls2recv[27];
    
    for ( t = 0; t < 27; t++ ) {
        sendBuf[t] = new double[domain2send[t]*PARTICLES_SIZE*sizeof(double)];
        recvBuf[t] = new double[EXPECTED_NUM_OF_PARTICLES*PARTICLES_SIZE*sizeof(double)];
        partcls2send[t] = 0;
        partcls2recv[t] = EXPECTED_NUM_OF_PARTICLES;
    }

    double* prtclPos;
    
    vector<int> removeFromLeaving;
    removeFromLeaving.reserve(leavingParticles.size());

    int posShift = phase == CORRECTOR ? 3 : 0;
    int outflowLeaving = 0;
    
    for ( int ptclNum = 0; ptclNum < leavingParticles.size(); ptclNum++ ){
        idx = leavingParticles[ptclNum];
        
        prtclPos = particles[idx]->getPosition();
        
        x0 = prtclPos[0+posShift];
        y0 = prtclPos[1+posShift];
        z0 = prtclPos[2+posShift];
        
        //check whom to send
        xloc = (x0 - domainShiftX)/domainXSize/dx;
        yloc = (y0 - domainShiftY)/domainYSize/dy;
        zloc = (z0 - domainShiftZ)/domainZSize/dz;
        
        switch (loader->dim) {
            case 1:
                a = (int) floor(xloc);
                b = 0;
                c = 0;
                break;
            case 2:
                a = (int) floor(xloc);
                b = (int) floor(yloc);
                c = 0;
                break;
            case 3:
                a = (int) floor(xloc);
                b = (int) floor(yloc);
                c = (int) floor(zloc);
                break;
                
            default:
                throw runtime_error("dimension is not 1/2/3");
        }
        
        t = (1+c)+3*((1+b)+3*(1+a));
 
#ifdef HEAVYLOG
        if ( t < 0 || t > 26 || std::isnan(x0) || std::isnan(y0) || std::isnan(z0) ){
            string pb ="[BoundaryManager] idx = "+to_string(idx)
            +"\n      x0 = "+to_string(x0)
            +"\n      y0 = "+to_string(y0)
            +"\n      z0 = "+to_string(z0)
            +"\n      xloc = "+to_string(xloc)
            +"\n      yloc = "+to_string(yloc)
            +"\n      zloc = "+to_string(zloc)
            +"\n      x1 = "+to_string(prtclPos[0+posShift])
            +"\n      y1 = "+to_string(prtclPos[1+posShift])
            +"\n      z1 = "+to_string(prtclPos[2+posShift])
            +"\n      a = "+to_string(a)
            +"\n      b = "+to_string(b)
            +"\n      c = "+to_string(c)
            +"\n      t = "+to_string(t);
            logger->writeMsg(pb.c_str(), CRITICAL);
        }
#endif
        
        if (  applyPeriodicBC(particles[idx], phase) == 1 ) {
            // need to send particle and need to remove from home domain
        }
        
//        if ( applyOutflowBC(particles[idx], phase) == 1 ) {
//            //no need to send particle but need to remove from home domain
//            // just keep it in leaving set and do not serialize
//            continue;
//        }
        if ( applyOutflowBC(t) == 1 ) {
            //no need to send particle but need to remove from home domain
            // just keep it in leaving set and do not serialize
            outflowLeaving++;
            continue;
        }

        
        // reflect BCs are out-of-date
//        if ( applyReflectBC(particles[idx], phase) == 1 ) {
//            //no need to send particle no need to remove from home domain
//            // remove from leaving set and do not serialize
//            removeFromLeaving.push_back(ptclNum);
//            continue;
//        }
        
        if( t != 13 ){
            particles[idx]->serialize(sendBuf[t], PARTICLES_SIZE*partcls2send[t]);
            partcls2send[t] += 1;
        } else {
            removeFromLeaving.push_back(ptclNum);
        }
    }
    
    for( int ptclNum = removeFromLeaving.size()-1; ptclNum >= 0 ; ptclNum-- ){
        leavingParticles.erase(leavingParticles.begin()+removeFromLeaving[ptclNum]);
    }
    
    string msd ="[BoundaryManager] total leaving = "+to_string(leavingParticles.size())
    + " / outflowLeaving = " + to_string(outflowLeaving);
    logger->writeMsg(msd.c_str(), DEBUG);
    
    MPI_Status st;
    int receivedTot;
    for ( t = 0; t < 27; t++ ){
            if(t != 13){
                
                MPI_Sendrecv(sendBuf[t], PARTICLES_SIZE*sizeof(double)*partcls2send[t],
                             MPI_DOUBLE, loader->neighbors2Send[t], t,
                             recvBuf[t], PARTICLES_SIZE*sizeof(double)*partcls2recv[t],
                             MPI_DOUBLE, loader->neighbors2Recv[t], t,
                             MPI_COMM_WORLD, &st);
                
                MPI_Get_count(&st, MPI_DOUBLE, &receivedTot);
                

                receivedTot = receivedTot/PARTICLES_SIZE/sizeof(double);
    
                for (int ptclNum = 0; ptclNum < receivedTot; ptclNum++){
                    particles2add.push_back(shared_ptr<Particle>(new Particle));
                    int idxOfAdded = particles2add.size()-1;
                    particles2add[idxOfAdded]->deserialize(recvBuf[t], PARTICLES_SIZE*ptclNum);
                }
                
            }
    }
    
    for ( t = 0; t < 27; t++ ) {
                delete [] sendBuf[t];
                delete [] recvBuf[t];
    }
    
    auto end_time = high_resolution_clock::now();
    string msg ="[BoundaryManager] apply BC duration = "
    +to_string(duration_cast<milliseconds>(end_time - start_time).count())+" ms";
    logger->writeMsg(msg.c_str(),  DEBUG);
}


int BoundaryManager::applyOutflowBC(Particle* particle, int phase){
    
    int remove = 0;
    double* prtclPos = particle->getPosition();
    int posShift = phase == CORRECTOR? 3 : 0;
    double L, r0;
    
    for (int i = 0; i < 3; i++) {
        if ( loader->partclBCtype[i] == OUTFLOW_BC ){
            L  = loader->boxSizes[i];
            r0 = prtclPos[i+posShift];
        
            if (r0 < 0.0){
                return 1;
            }
            if (r0 > L){
                return 1;
            }
        }
    }
    return remove;
}

int BoundaryManager::applyOutflowBC(int sendTo){
    
    int remove = 0;
    
    if ( loader->partclBCtype[0] == OUTFLOW_BC ){
        if ( sendTo == NEIGHBOR_LEFT && loader->neighbors2Send[NEIGHBOR_LEFT] == MPI_PROC_NULL ){
            remove = 1;
        } else if ( sendTo == NEIGHBOR_RIGHT && loader->neighbors2Send[NEIGHBOR_RIGHT] == MPI_PROC_NULL){
            remove = 1;
        }
    }
    
    if ( loader->partclBCtype[1] == OUTFLOW_BC ){
        if ( sendTo == NEIGHBOR_BOTTOM && loader->neighbors2Send[NEIGHBOR_BOTTOM] == MPI_PROC_NULL){
            remove = 1;
        } else if ( sendTo == NEIGHBOR_TOP && loader->neighbors2Send[NEIGHBOR_TOP] == MPI_PROC_NULL){
            remove = 1;
        }
    }
    
    if ( loader->partclBCtype[2] == OUTFLOW_BC ){
        if ( sendTo == NEIGHBOR_BACK && loader->neighbors2Send[NEIGHBOR_BACK] == MPI_PROC_NULL){
            remove = 1;
        } else if ( sendTo == NEIGHBOR_FRONT && loader->neighbors2Send[NEIGHBOR_FRONT] == MPI_PROC_NULL){
            remove = 1;
        }
    }
    
    return remove;
}





//shift coordinates
int BoundaryManager::applyPeriodicBC(Particle* particle, int phase){
    
    double* prtclPos = particle->getPosition();
    
    int posShift = phase == CORRECTOR? 3 : 0;
    
    double L, r0, rnew;
    int period = 0;
    for (int i = 0; i < 3; i++) {
        L  = loader->boxSizes[i];
        r0 = prtclPos[i+posShift];
        
        if ( loader->partclBCtype[i] == PERIODIC_BC ){

            if (r0 < 0.0){
                rnew = (r0+L) > L ? L-EPS4 : r0+L;
                particle->setPosition(i,   rnew);
                particle->setPosition(i+3, rnew);
                period = 1;
            }
            if (r0 >= L){
                rnew = (r0-L) < 0 ? EPS4 : r0-L;
                particle->setPosition(i,   rnew);
                particle->setPosition(i+3, rnew);
                period = 1;
            }
        }
    }
    return period;
}





//reflect coordinates
int BoundaryManager::applyReflectBC(Particle* particle, int phase){
    
    double* prtclPos = particle->getPosition();
    double* prtclVel = particle->getVelocity();
    
    int shift = phase == CORRECTOR? 3 : 0;
    
    double L, r0, V0;
    int revert = 0;
    for ( int i = 0; i < 3; i++) {
        if ( loader->partclBCtype[i] == REFLECT_BC ){
            
            L  = loader->boxSizes[i];
            r0 = prtclPos[i+shift];
            V0 = prtclVel[i+shift];
            
            if ( r0 < 0.0 ) {
                particle->setPosition(i+3, -r0);
                particle->setVelocity(i+3, -V0);
                particle->setPosition(i, -r0);
                particle->setVelocity(i, -V0);
                revert = 1;
            }
            
            if ( r0 >= L ) {
                particle->setPosition(i+3, -r0 + 2*L);
                particle->setVelocity(i+3, -V0);
                particle->setPosition(i, -r0 + 2*L);
                particle->setVelocity(i, -V0);
                revert = 1;
            }
           
        }
    }
    
    return revert;
}


vector<int> BoundaryManager::getLeavingParticlesIdxs(){
    return leavingParticles;
}

void BoundaryManager::reset(){
    logger->writeMsg("[BoundaryManager] reset ..", DEBUG);
    leavingParticles.clear();
    
    for ( int t = 0; t < 27; t++ ) {
        domain2send[t] = 0;
    }
}

void BoundaryManager::storeParticle(int idx, double pos[3]){
    double domainShiftX = loader->boxCoordinates[0][0];
    double domainShiftY = loader->boxCoordinates[1][0];
    double domainShiftZ = loader->boxCoordinates[2][0];
    double domainXSize = loader->resolution[0];
    double domainYSize = loader->resolution[1];
    double domainZSize = loader->resolution[2];
    double dx = loader->spatialSteps[0];
    double dy = loader->spatialSteps[1];
    double dz = loader->spatialSteps[2];
    
    double x, y, z;
    int a, b, c, domain;
    
    x = (pos[0] - domainShiftX)/domainXSize/dx;
    y = (pos[1] - domainShiftY)/domainYSize/dy;
    z = (pos[2] - domainShiftZ)/domainZSize/dz;
    
  
    switch (loader->dim) {
        case 1:
            a = (int)floor(x);
            b = 0;
            c = 0;
            break;
        case 2:
            a = (int)floor(x);
            b = (int)floor(y);
            c = 0;
            break;
        case 3:
            a = (int)floor(x);
            b = (int)floor(y);
            c = (int)floor(z);
            break;
            
        default:
            throw runtime_error("dimension is not 1/2/3");
    }

    
    domain = (1+c)+3*((1+b)+3*(1+a));

    leavingParticles.push_back(idx);

    
    if ( domain < 0 || domain > 26){
        #ifdef HEAVYLOG
        string pb ="[BoundaryManager] problem to store idx = "+to_string(idx)
        +"\n      x0 = "+to_string(x)
        +"\n      y0 = "+to_string(y)
        +"\n      z0 = "+to_string(z)
        +"\n      pos[0] = "+to_string(pos[0])
        +"\n      pos[1] = "+to_string(pos[1])
        +"\n      pos[2] = "+to_string(pos[2])
        +"\n      a = "+to_string(a)
        +"\n      b = "+to_string(b)
        +"\n      c = "+to_string(c)
        +"\n      domain = "+to_string(domain);
        logger->writeMsg(pb.c_str(), CRITICAL);
        #endif
    }else{
        // need to know before sending for memory preallocation
        domain2send[domain] += 1;
    }
    

}












