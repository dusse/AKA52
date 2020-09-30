#include "GridManager.hpp"

using namespace std;
using namespace chrono;

GridManager::GridManager(shared_ptr<Loader> loader){
    this->loader=loader;
    logger.reset(new Logger());
    initialize();
    string msg ="[GridManager] init...OK {Node number:"+to_string(totalNodeNumber)+"}";
    logger->writeMsg(msg.c_str(), DEBUG);
}


GridManager::~GridManager(){
    delete [] neibors4G2spatialDerX;
    delete [] neibors4G2spatialDerY;
    delete [] neibors4G2spatialDerZ;
    delete [] negbors4PresX;
    delete [] negbors4PresY;
    delete [] negbors4PresZ;
    delete [] negbors4LaplacX;
    delete [] negbors4LaplacY;
    delete [] negbors4LaplacZ;
    delete [] neibors4G1spatialDerX;
    delete [] neibors4G1spatialDerY;
    delete [] neibors4G1spatialDerZ;
    delete [] neighbourhood;
    delete [] nodesG2vars;
    delete [] nodesG1vars;
    
    for (int t = 0; t < 27; t++) {
        delete [] sendIdx4Gath[t];
        delete [] recvIdx4Gath[t];
        delete [] sendIdxOnG4[t];
        delete [] recvIdxOnG4[t];
        delete [] sendIdx[t];
        delete [] recvIdx[t];
    }
}

const int* GridManager::getNeighbourhoodOnG1(){
    return neighbourhood;
}

void GridManager::initialize(){
    logger->writeMsg("[GridManager] initialize() ...", DEBUG);
    
    int xRes = loader->resolution[0];
    int yRes = loader->resolution[1];
    int zRes = loader->resolution[2];
    
    this->totalNodeNumber = xRes*yRes*zRes;
    
    this->G1nodesNumber = (xRes+1)*(yRes+1)*(zRes+1);
    this->G2nodesNumber = (xRes+2)*(yRes+2)*(zRes+2);
    
     neibors4G2spatialDerX = new int[G2nodesNumber*8*sizeof(int)];
     neibors4G2spatialDerY = new int[G2nodesNumber*8*sizeof(int)];
     neibors4G2spatialDerZ = new int[G2nodesNumber*8*sizeof(int)];
    
     negbors4PresX = new int[G2nodesNumber*18*sizeof(int)];
     negbors4PresY = new int[G2nodesNumber*18*sizeof(int)];
     negbors4PresZ = new int[G2nodesNumber*18*sizeof(int)];
    
     negbors4LaplacX = new int[G2nodesNumber*4*sizeof(int)];
     negbors4LaplacY = new int[G2nodesNumber*4*sizeof(int)];
     negbors4LaplacZ = new int[G2nodesNumber*4*sizeof(int)];
    
     neibors4G1spatialDerX = new int[G1nodesNumber*8*sizeof(int)];
     neibors4G1spatialDerY = new int[G1nodesNumber*8*sizeof(int)];
     neibors4G1spatialDerZ = new int[G1nodesNumber*8*sizeof(int)];

    neighbourhood = new int[G1nodesNumber*8*sizeof(int)];

    int numOfSpecies = loader->getNumberOfSpecies();
    
    nodesG2vars = new VectorVar*[G2nodesNumber*(NUM_OF_MAIN_G2VARS+2*numOfSpecies)];
    nodesG1vars = new VectorVar*[G1nodesNumber*2];
    
    initG1Nodes();
    initG2Nodes();

    sendRecvIndecis4MPI();
    sendRecvIndecis4MPIext();
    sendRecvIndecis4MPIonG4();
    
    initBoundaryIndecies();
    
    getVectorVariablesForAllNodes();
    logger->writeMsg("[GridManager] initialize() ...OK", DEBUG);
}

const int* GridManager::getNeibors4G2spatialDerX(){
    return neibors4G2spatialDerX;
}

const int* GridManager::getNeibors4G2spatialDerY(){
    return neibors4G2spatialDerY;
}

const int* GridManager::getNeibors4G2spatialDerZ(){
    return neibors4G2spatialDerZ;
}


const int* GridManager::getNeibors4G1spatialDerX(){
    return neibors4G1spatialDerX;
}

const int* GridManager::getNeibors4G1spatialDerY(){
    return neibors4G1spatialDerY;
}

const int* GridManager::getNeibors4G1spatialDerZ(){
    return neibors4G1spatialDerZ;
}


int GridManager::getVarsNumOnG1(){
    return totVarsOnG1;
}


/*                        E(i,j+1)           |              E(i+1,j+1)
 *                           x-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-x
 *                           .               |               .
 *                           |               |   CELL i,j    |
 *                           .               |               .
 *                         dy/2              |   O(x,y)      |
 *                           .               |     ->        .
 *                           |               |i=int(x/dx+0.5)|
 *                           .               |j=int(y/dy+0.5).
 *                           |               |               |
 *            .B(i-1,j)______________________.B(i,j)____________________.B(i+1,j)
 *                           |               |               |
 *                           .               |               .
 *                           |               |               |
 *                         dy/2              |               .
 *                           |               |               |
 *                           .               |               .
 *                           |               |               |
 *                           . E(i,j)        |               .
 *                           x-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-x E(i+1,j)
 *                                  dx/2     |      dx/2
 */

#define NEIGHBOR_LEFT   4
#define NEIGHBOR_RIGHT  22
#define NEIGHBOR_BOTTOM 10
#define NEIGHBOR_TOP    16
#define NEIGHBOR_BACK   12
#define NEIGHBOR_FRONT  14

void GridManager::initBoundaryIndecies(){
    int xRes = loader->resolution[0];
    int yRes = loader->resolution[1];
    int zRes = loader->resolution[2];
    int xResG2 = xRes+2, yResG2 = yRes+2, zResG2 = zRes+2;
    int i,j,k, idx2set, idx2use;

    if(loader->BCtype[0] == IDEAL && xRes != 1){
        
        if (loader->neighbors2Send[NEIGHBOR_LEFT] == MPI_PROC_NULL){
            for ( j=0; j<yResG2; j++){
                for ( k=0; k<zResG2; k++){
                    idx2use = IDX(1,j,k, xResG2,yResG2,zResG2);
                    idx2set = IDX(0,j,k, xResG2,yResG2,zResG2);
                        
                    idxs4BoundaryX2fill[idx2use] = idx2set;
                }
            }
        }
        
        if (loader->neighbors2Send[NEIGHBOR_RIGHT] == MPI_PROC_NULL){
            for ( j=0; j<yResG2; j++){
                for ( k=0; k<zResG2; k++){
                    idx2use = IDX(xResG2-2,j,k, xResG2,yResG2,zResG2);
                    idx2set = IDX(xResG2-1,j,k, xResG2,yResG2,zResG2);
                    
                    idxs4BoundaryX2fill[idx2use] = idx2set;
                }
            }
        }
    }
    
    if(loader->BCtype[1] == IDEAL && yRes != 1){
        
        if (loader->neighbors2Send[NEIGHBOR_BOTTOM] == MPI_PROC_NULL){
            for ( i=0; i<xResG2; i++){
                for ( k=0; k<zResG2; k++){
                    idx2use = IDX(i,1,k, xResG2,yResG2,zResG2);
                    idx2set = IDX(i,0,k, xResG2,yResG2,zResG2);
                    
                    idxs4BoundaryY2fill[idx2use] = idx2set;
                }
            }
        }
        if (loader->neighbors2Send[NEIGHBOR_TOP] == MPI_PROC_NULL){
            for ( i=0; i<xResG2; i++){
                for ( k=0; k<zResG2; k++){
                    idx2use = IDX(i,yResG2-2,k, xResG2,yResG2,zResG2);
                    idx2set = IDX(i,yResG2-1,k, xResG2,yResG2,zResG2);
                    
                    idxs4BoundaryY2fill[idx2use] = idx2set;
                }
            }
        }
    }
    
    if(loader->BCtype[2] == IDEAL && zRes != 1){
        
        if (loader->neighbors2Send[NEIGHBOR_BACK] == MPI_PROC_NULL){
            for ( i=0; i<xResG2; i++){
                for ( j=0; j<yResG2; j++){
                    idx2use = IDX(i,j,1, xResG2,yResG2,zResG2);
                    idx2set = IDX(i,j,0, xResG2,yResG2,zResG2);
                    
                    idxs4BoundaryZ2fill[idx2use] = idx2set;
                }
            }
        }
        if (loader->neighbors2Send[NEIGHBOR_FRONT] == MPI_PROC_NULL){
            for ( i=0; i<xResG2; i++){
                for ( j=0; j<yResG2; j++){
                    idx2use = IDX(i,j,zResG2-2, xResG2,yResG2,zResG2);
                    idx2set = IDX(i,j,zResG2-1, xResG2,yResG2,zResG2);
                    
                    idxs4BoundaryZ2fill[idx2use] = idx2set;
                }
            }
        }
    }
    
    string msg ="[GridManager] initialized idxs4BoundaryX2fill.size() = "
    +to_string(idxs4BoundaryX2fill.size())
    +"\n                       idxs4BoundaryY2fill.size() = "
    +to_string(idxs4BoundaryY2fill.size())
    +"\n                       idxs4BoundaryZ2fill.size() = "
    +to_string(idxs4BoundaryZ2fill.size());
    logger->writeMsg(msg.c_str(), DEBUG);
}

//G1 grid describes vecticies of the cells
void GridManager::initG1Nodes(){
    auto start_time = high_resolution_clock::now();
    logger->writeMsg("[GridManager] init Node on G1 start ", DEBUG);
    int xRes = loader->resolution[0];
    int yRes = loader->resolution[1];
    int zRes = loader->resolution[2];
    int xResG1 = xRes+1, yResG1 = yRes+1, zResG1 = zRes+1;
    int xResG2 = xRes+2, yResG2 = yRes+2, zResG2 = zRes+2;
    
    int   dx[4][2][3] = {{{1,0,0}, {0,0,0}}, {{1,1,0}, {0,1,0}},
                         {{1,0,1}, {0,0,1}}, {{1,1,1}, {0,1,1}}};
    int   dy[4][2][3] = {{{0,1,0}, {0,0,0}}, {{1,1,0}, {1,0,0}},
                         {{0,1,1}, {0,0,1}}, {{1,1,1}, {1,0,1}}};
    int   dz[4][2][3] = {{{0,0,1}, {0,0,0}}, {{1,0,1}, {1,0,0}},
                         {{0,1,1}, {0,1,0}}, {{1,1,1}, {1,1,0}}};
    
    int neibIDXs[8][3] = {{0,0,0}, {0,0,1}, {0,1,0}, {1,0,0},
                          {0,1,1}, {1,1,0}, {1,0,1}, {1,1,1}};
    int i,j,k,idx;
    int idx_count = 0;
    for ( i=0; i<xResG1; i++){
        for ( j=0; j<yResG1; j++){
            for ( k=0; k<zResG1; k++){
                idx   = IDX(i,j,k,xResG1,yResG1,zResG1);
                
                for(int comp = 0; comp < 8; comp++){
                    neibors4G2spatialDerX[8*idx+comp] = 0;
                    neibors4G2spatialDerY[8*idx+comp] = 0;
                    neibors4G2spatialDerZ[8*idx+comp] = 0;
                    neighbourhood[8*idx+comp] = 0;
                }

                
                vector<int> neiborsX = {0, 0, 0, 0, 0, 0, 0, 0};
                vector<int> neiborsY = {0, 0, 0, 0, 0, 0, 0, 0};
                vector<int> neiborsZ = {0, 0, 0, 0, 0, 0, 0, 0};
                    
                for(int pairNum = 0; pairNum<4; pairNum++){
                        if(xRes != 1){
                            neiborsX[2*pairNum+0] = IDX(i+dx[pairNum][0][0],
                                                        j+dx[pairNum][0][1],
                                                        k+dx[pairNum][0][2],
                                                        xResG2,yResG2,zResG2);
                            
                            neiborsX[2*pairNum+1] = IDX(i+dx[pairNum][1][0],
                                                        j+dx[pairNum][1][1],
                                                        k+dx[pairNum][1][2],
                                                        xResG2,yResG2,zResG2);
                        }
                        if(yRes != 1){
                            neiborsY[2*pairNum+0] = IDX(i+dy[pairNum][0][0],
                                                        j+dy[pairNum][0][1],
                                                        k+dy[pairNum][0][2],
                                                        xResG2,yResG2,zResG2);
                            
                            neiborsY[2*pairNum+1] = IDX(i+dy[pairNum][1][0],
                                                        j+dy[pairNum][1][1],
                                                        k+dy[pairNum][1][2],
                                                        xResG2,yResG2,zResG2);
                        }
                        
                        if(zRes != 1){
                            neiborsZ[2*pairNum+0] = IDX(i+dz[pairNum][0][0],
                                                        j+dz[pairNum][0][1],
                                                        k+dz[pairNum][0][2],
                                                        xResG2,yResG2,zResG2);
                            
                            neiborsZ[2*pairNum+1] = IDX(i+dz[pairNum][1][0],
                                                        j+dz[pairNum][1][1],
                                                        k+dz[pairNum][1][2],
                                                        xResG2,yResG2,zResG2);
                        }
                        
                }
                    
                for(int comp = 0; comp < 8; comp++){
                    neibors4G2spatialDerX[8*idx+comp] = neiborsX[comp];
                    neibors4G2spatialDerY[8*idx+comp] = neiborsY[comp];
                    neibors4G2spatialDerZ[8*idx+comp] = neiborsZ[comp];
                }
                    
                if(i != 0 && j != 0 && k != 0){
                    for(int comp = 0; comp < 8; comp++){
                        neiborsX[comp] = 0;
                        neiborsY[comp] = 0;
                        neiborsZ[comp] = 0;
                    }
                    
                    for(int pairNum = 0; pairNum<4; pairNum++){
                            if(xRes != 1){
                                neiborsX[2*pairNum+0] = IDX(i-dx[pairNum][0][0],
                                                            j-dx[pairNum][0][1],
                                                            k-dx[pairNum][0][2],
                                                            xResG1,yResG1,zResG1);
                                
                                neiborsX[2*pairNum+1] = IDX(i-dx[pairNum][1][0],
                                                            j-dx[pairNum][1][1],
                                                            k-dx[pairNum][1][2],
                                                            xResG1,yResG1,zResG1);
                            }
                            if(yRes != 1){
                                neiborsY[2*pairNum+0] = IDX(i-dy[pairNum][0][0],
                                                            j-dy[pairNum][0][1],
                                                            k-dy[pairNum][0][2],
                                                            xResG1,yResG1,zResG1);
                                
                                neiborsY[2*pairNum+1] = IDX(i-dy[pairNum][1][0],
                                                            j-dy[pairNum][1][1],
                                                            k-dy[pairNum][1][2],
                                                            xResG1,yResG1,zResG1);
                            }
                            if(zRes != 1){
                                neiborsZ[2*pairNum+0] = IDX(i-dz[pairNum][0][0],
                                                            j-dz[pairNum][0][1],
                                                            k-dz[pairNum][0][2],
                                                            xResG1,yResG1,zResG1);
                                
                                neiborsZ[2*pairNum+1] = IDX(i-dz[pairNum][1][0],
                                                            j-dz[pairNum][1][1],
                                                            k-dz[pairNum][1][2],
                                                            xResG1,yResG1,zResG1);
                            }
                    }
                        
                    for(int comp = 0; comp < 8; comp++){
                        neibors4G1spatialDerX[8*idx+comp] = neiborsX[comp];
                        neibors4G1spatialDerY[8*idx+comp] = neiborsY[comp];
                        neibors4G1spatialDerZ[8*idx+comp] = neiborsZ[comp];
                    }
                    
                    idx_count++;
                    
                    int a = xRes == 1 ? 0 : 1;
                    int b = yRes == 1 ? 0 : 1;
                    int c = zRes == 1 ? 0 : 1;
                    
                    for(int neiNum = 0; neiNum<8; neiNum++){
                        neighbourhood[8*idx+neiNum] = IDX(i-a*neibIDXs[neiNum][0],
                                                          j-b*neibIDXs[neiNum][1],
                                                          k-c*neibIDXs[neiNum][2],
                                                          xResG1,yResG1,zResG1);
                    }
                }
                
                VectorVar* mag     = new VectorVar(MAGNETIC,     {0.0, 0.0, 0.0});
                VectorVar* mag_aux = new VectorVar(MAGNETIC_AUX, {0.0, 0.0, 0.0});
                nodesG1vars[G1nodesNumber*MAGNETIC     + idx] = mag;
                nodesG1vars[G1nodesNumber*MAGNETIC_AUX + idx] = mag_aux;
                totVarsOnG1 = SIZEG1;
            }
        }
    }
    
   
    
    auto end_time = high_resolution_clock::now();
    string msg ="[GridManager] init Node G1 OK..duration = "
                +to_string(duration_cast<milliseconds>(end_time - start_time).count())+" ms";
    logger->writeMsg(msg.c_str(), DEBUG);
}

const int* GridManager::getNghbd4DivPeInXonG2(){
    return negbors4PresX;
}

const int* GridManager::getNghbd4DivPeInYonG2(){
    return negbors4PresY;
}

const int* GridManager::getNghbd4DivPeInZonG2(){
    return negbors4PresZ;
}

const int* GridManager::getNghbd4LaplacInXonG2(){
    return negbors4LaplacX;
}

const int* GridManager::getNghbd4LaplacInYonG2(){
    return negbors4LaplacY;
}

const int* GridManager::getNghbd4LaplacInZonG2(){
    return negbors4LaplacZ;
}


int GridManager::getVarsNumOnG2(){
    return totVarsOnG2;
}


//G2 grid describes extended grid G1
void GridManager::initG2Nodes(){
    auto start_time = high_resolution_clock::now();
    logger->writeMsg("[GridManager] init Node start ", DEBUG);
    int xRes = loader->resolution[0], yRes = loader->resolution[1], zRes = loader->resolution[2];
    int xResG2 = xRes+2, yResG2 = yRes+2, zResG2 = zRes+2;
//                                    +          -
    int PAIRS4divPinX[9][2][3] = { {{ 1, 1,-1},{-1, 1,-1}},
                                   {{ 1,-1,-1},{-1,-1,-1}}, // weight  = 0.125
                                   {{ 1, 1, 1},{-1, 1, 1}},
                                   {{ 1,-1, 1},{-1,-1, 1}},
        
                                   {{ 1, 0,-1},{-1, 0,-1}},
                                   {{ 1, 1, 0},{-1, 1, 0}}, // weight  = 0.25
                                   {{ 1,-1, 0},{-1,-1, 0}},
                                   {{ 1, 0, 1},{-1, 0, 1}},
        
                                   {{ 1, 0, 0},{-1, 0, 0}}};// weight  = 0.5
//                                       +          -
    int PAIRS4divPinY[9][2][3] = { {{ 1, 1,-1},{ 1,-1,-1}},
                                   {{-1, 1,-1},{-1,-1,-1}},
                                   {{ 1, 1, 1},{ 1,-1, 1}},
                                   {{-1, 1, 1},{-1,-1, 1}},
        
                                   {{ 0, 1,-1},{ 0,-1,-1}},
                                   {{ 1, 1, 0},{ 1,-1, 0}},
                                   {{-1, 1, 0},{-1,-1, 0}},
                                   {{ 0, 1, 1},{ 0,-1, 1}},
        
                                   {{ 0, 1, 0},{ 0,-1, 0}}};
//                                          +          -
    int PAIRS4divPinZ[9][2][3] = { {{ 1,-1, 1},{ 1,-1,-1}},
                                   {{-1,-1, 1},{-1,-1,-1}},
                                   {{ 1, 1, 1},{ 1, 1,-1}},
                                   {{-1, 1, 1},{-1, 1,-1}},
       
                                   {{ 0,-1, 1},{ 0,-1,-1}},
                                   {{ 1, 0, 1},{ 1, 0,-1}},
                                   {{-1, 0, 1},{-1, 0,-1}},
                                   {{ 0, 1, 1},{ 0, 1,-1}},
        
                                   {{ 0, 0, 1},{ 0, 0,-1}}};
    
    int PAIRS4LapJinX[2][2][3] = {{{1,0,0}, {0,0,0}}, {{-1,0,0}, {0,0,0}}};
    int PAIRS4LapJinY[2][2][3] = {{{0,1,0}, {0,0,0}}, {{0,-1,0}, {0,0,0}}};
    int PAIRS4LapJinZ[2][2][3] = {{{0,0,1}, {0,0,0}}, {{0,0,-1}, {0,0,0}}};

    int numOfSpecies = loader->getNumberOfSpecies();

    int i,j,k,idx, idxG1;
    for ( i=0; i<xResG2; i++){
        for ( j=0; j<yResG2; j++){
            for ( k=0; k<zResG2; k++){
                idx = IDX(i,j,k,xResG2,yResG2,zResG2);
                
                vector<int> neibors4divPX(18);
                vector<int> neibors4divPY(18);
                vector<int> neibors4divPZ(18);
                
                for(int comp = 0; comp < 18; comp++){
                     negbors4PresX[18*idx+comp] = 0;
                     negbors4PresY[18*idx+comp] = 0;
                     negbors4PresZ[18*idx+comp] = 0;
                }
                
                vector<int> neibors4LapX(4);
                vector<int> neibors4LapY(4);
                vector<int> neibors4LapZ(4);
                
                for(int comp = 0; comp < 4; comp++){
                    negbors4LaplacX[4*idx+comp] = 0;
                    negbors4LaplacY[4*idx+comp] = 0;
                    negbors4LaplacZ[4*idx+comp] = 0;
                }
                
                if(i != 0 && j != 0 && k != 0 && i != xResG2-1
                          && j != yResG2-1 && k != zResG2-1){
                    
                    for(int pairNum = 0; pairNum<9; pairNum++){
                        if(xRes != 1){
                            neibors4divPX[2*pairNum+0] = IDX(i+PAIRS4divPinX[pairNum][0][0],
                                                             j+PAIRS4divPinX[pairNum][0][1],
                                                             k+PAIRS4divPinX[pairNum][0][2],
                                                             xResG2,yResG2,zResG2);
                            
                            neibors4divPX[2*pairNum+1] = IDX(i+PAIRS4divPinX[pairNum][1][0],
                                                             j+PAIRS4divPinX[pairNum][1][1],
                                                             k+PAIRS4divPinX[pairNum][1][2],
                                                             xResG2,yResG2,zResG2);
                        }
                        if(yRes != 1){
                            neibors4divPY[2*pairNum+0] = IDX(i+PAIRS4divPinY[pairNum][0][0],
                                                             j+PAIRS4divPinY[pairNum][0][1],
                                                             k+PAIRS4divPinY[pairNum][0][2],
                                                             xResG2,yResG2,zResG2);
                            
                            neibors4divPY[2*pairNum+1] = IDX(i+PAIRS4divPinY[pairNum][1][0],
                                                             j+PAIRS4divPinY[pairNum][1][1],
                                                             k+PAIRS4divPinY[pairNum][1][2],
                                                             xResG2,yResG2,zResG2);
                        }
                        if(zRes != 1){
                            neibors4divPZ[2*pairNum+0] = IDX(i+PAIRS4divPinZ[pairNum][0][0],
                                                             j+PAIRS4divPinZ[pairNum][0][1],
                                                             k+PAIRS4divPinZ[pairNum][0][2],
                                                             xResG2,yResG2,zResG2);
                        
                            neibors4divPZ[2*pairNum+1] = IDX(i+PAIRS4divPinZ[pairNum][1][0],
                                                             j+PAIRS4divPinZ[pairNum][1][1],
                                                             k+PAIRS4divPinZ[pairNum][1][2],
                                                             xResG2,yResG2,zResG2);
                        }
                        
                    }
                    
                    for(int comp = 0; comp < 18; comp++){
                        negbors4PresX[18*idx+comp] = neibors4divPX[comp];
                        negbors4PresY[18*idx+comp] = neibors4divPY[comp];
                        negbors4PresZ[18*idx+comp] = neibors4divPZ[comp];
                    }
                    
                    for(int pairNum = 0; pairNum<2; pairNum++){
                        if(xRes != 1){
                            neibors4LapX[2*pairNum+0] = IDX(i+PAIRS4LapJinX[pairNum][0][0],
                                                            j+PAIRS4LapJinX[pairNum][0][1],
                                                            k+PAIRS4LapJinX[pairNum][0][2],
                                                            xResG2,yResG2,zResG2);
                            
                            neibors4LapX[2*pairNum+1] = IDX(i+PAIRS4LapJinX[pairNum][1][0],
                                                            j+PAIRS4LapJinX[pairNum][1][1],
                                                            k+PAIRS4LapJinX[pairNum][1][2],
                                                            xResG2,yResG2,zResG2);
                        }
                        if(yRes != 1){
                            neibors4LapY[2*pairNum+0] = IDX(i+PAIRS4LapJinY[pairNum][0][0],
                                                            j+PAIRS4LapJinY[pairNum][0][1],
                                                            k+PAIRS4LapJinY[pairNum][0][2],
                                                            xResG2,yResG2,zResG2);
                            
                            neibors4LapY[2*pairNum+1] = IDX(i+PAIRS4LapJinY[pairNum][1][0],
                                                            j+PAIRS4LapJinY[pairNum][1][1],
                                                            k+PAIRS4LapJinY[pairNum][1][2],
                                                            xResG2,yResG2,zResG2);
                        }
                        if(zRes != 1){
                            neibors4LapZ[2*pairNum+0] = IDX(i+PAIRS4LapJinZ[pairNum][0][0],
                                                            j+PAIRS4LapJinZ[pairNum][0][1],
                                                            k+PAIRS4LapJinZ[pairNum][0][2],
                                                            xResG2,yResG2,zResG2);
                            
                            neibors4LapZ[2*pairNum+1] = IDX(i+PAIRS4LapJinZ[pairNum][1][0],
                                                            j+PAIRS4LapJinZ[pairNum][1][1],
                                                            k+PAIRS4LapJinZ[pairNum][1][2],
                                                            xResG2,yResG2,zResG2);
                        }
                    }
                    
                    for(int comp = 0; comp < 4; comp++){
                        negbors4LaplacX[4*idx+comp] = neibors4LapX[comp];
                        negbors4LaplacY[4*idx+comp] = neibors4LapY[comp];
                        negbors4LaplacZ[4*idx+comp] = neibors4LapZ[comp];
                    }

                }
                
                
                VectorVar* ele         = new VectorVar(ELECTRIC,     {0.0, 0.0, 0.0});
                VectorVar* ele_aux     = new VectorVar(ELECTRIC_AUX, {0.0, 0.0, 0.0});
                VectorVar* current     = new VectorVar(CURRENT,      {0.0, 0.0, 0.0});
                VectorVar* current_aux = new VectorVar(CURRENT_AUX,  {0.0, 0.0, 0.0});
                VectorVar* veloion     = new VectorVar(VELOCION,     {0.0, 0.0, 0.0});
                VectorVar* densele     = new VectorVar(DENSELEC,     {0.0});
                VectorVar* veloele     = new VectorVar(VELOCELE,     {0.0, 0.0, 0.0});
                VectorVar* presure     = new VectorVar(PRESSURE,     {0.0, 0.0, 0.0, 0.0, 0.0, 0.0});
                VectorVar* presure_aux = new VectorVar(PRESSURE_AUX, {0.0, 0.0, 0.0, 0.0, 0.0, 0.0});
                VectorVar* presure_smo = new VectorVar(PRESSURE_SMO, {0.0, 0.0, 0.0, 0.0, 0.0, 0.0});
                VectorVar* pdriver     = new VectorVar(DRIVER,       {0.0, 0.0, 0.0, 0.0, 0.0, 0.0});
                VectorVar* pdriver_aux = new VectorVar(DRIVER_AUX,   {0.0, 0.0, 0.0, 0.0, 0.0, 0.0});
                
                totVarsOnG2 = NUM_OF_MAIN_G2VARS+2*numOfSpecies;
                
                nodesG2vars[G2nodesNumber*ELECTRIC+    idx] = ele;
                nodesG2vars[G2nodesNumber*ELECTRIC_AUX+idx] = ele_aux;
                nodesG2vars[G2nodesNumber*CURRENT     +idx] = current;
                nodesG2vars[G2nodesNumber*CURRENT_AUX +idx] = current_aux;
                nodesG2vars[G2nodesNumber*VELOCION    +idx] = veloion;
                nodesG2vars[G2nodesNumber*DENSELEC    +idx] = densele;
                nodesG2vars[G2nodesNumber*VELOCELE    +idx] = veloele;
                nodesG2vars[G2nodesNumber*PRESSURE    +idx] = presure;
                nodesG2vars[G2nodesNumber*PRESSURE_AUX+idx] = presure_aux;
                nodesG2vars[G2nodesNumber*PRESSURE_SMO+idx] = presure_smo;
                nodesG2vars[G2nodesNumber*DRIVER      +idx] = pdriver;
                nodesG2vars[G2nodesNumber*DRIVER_AUX  +idx] = pdriver_aux;
                
                SHIFT_MAIN_DENS = NUM_OF_MAIN_G2VARS;
                SHIFT_MAIN_DENS_AUX  = SHIFT_MAIN_DENS+numOfSpecies;
                
                for(int spn = 0; spn<numOfSpecies;spn++){
                    VectorVar* densvel = new VectorVar(DENS_VEL(spn), {0.0, 0.0, 0.0, 0.0});
                    VectorVar* densaux = new VectorVar(DENS_AUX(spn), {0.0});
                    
                    nodesG2vars[G2nodesNumber*(SHIFT_MAIN_DENS    +spn)+idx] = densvel;
                    nodesG2vars[G2nodesNumber*(SHIFT_MAIN_DENS_AUX+spn)+idx] = densaux;
                }
                
            }
        }
    }
    auto end_time = high_resolution_clock::now();
    string msg ="[GridManager] init Node G2 duration = "
                +to_string(duration_cast<milliseconds>(end_time - start_time).count())+" ms";
    logger->writeMsg(msg.c_str(), DEBUG);
}

int GridManager::DENS_AUX(int sp){
    return SHIFT_MAIN_DENS_AUX+sp;
}

int GridManager::DENS_VEL(int sp){
    return SHIFT_MAIN_DENS+sp;
}


void GridManager::applyBC(int varName){
    //nothing for periodic BC, maps are empty
    auto start_time = high_resolution_clock::now();
    
    int varDim = nodesG2vars[G2nodesNumber*varName]->getSize();
    int varShift = G2nodesNumber*varName;
    int dim;
    
    const double* vectorVar;
    
    map<int, int>  allBoundaryCells;
    allBoundaryCells.insert(idxs4BoundaryX2fill.begin(), idxs4BoundaryX2fill.end());
    allBoundaryCells.insert(idxs4BoundaryY2fill.begin(), idxs4BoundaryY2fill.end());
    allBoundaryCells.insert(idxs4BoundaryZ2fill.begin(), idxs4BoundaryZ2fill.end());
    
    set<int> specialTreat;
    int numOfSpecies = loader->getNumberOfSpecies();
    for(int i=0;i<numOfSpecies;i++){
        specialTreat.insert(DENS_VEL(i));
    }
    
    if ( specialTreat.count(varName) ){
        for ( const auto &keyval : allBoundaryCells ) {
            vectorVar = nodesG2vars[varShift+keyval.first]->getValue();
            for ( dim=0; dim<varDim; dim++) {
                double oldVal2use = vectorVar[dim];
                double oldVal2set = nodesG2vars[varShift+keyval.second]->getValue()[dim];
                nodesG2vars[varShift+keyval.second]->addValue( dim, oldVal2use);
                nodesG2vars[varShift+keyval.first]->addValue( dim, oldVal2set);
            }
        }
        auto end_time = high_resolution_clock::now();
        string msg ="[GridManager] apply BC for var = "+to_string(varName) +"  duration = "
        +to_string(duration_cast<milliseconds>(end_time - start_time).count())+" ms";
        logger->writeMsg(msg.c_str(), DEBUG);
        
        return;
    }


    set<int> simpleRevertVariables = {CURRENT, CURRENT_AUX, VELOCION};
    
    if ( simpleRevertVariables.count(varName) ){
        
        for ( const auto &keyval : idxs4BoundaryX2fill ) {
            vectorVar = nodesG2vars[varShift+keyval.first]->getValue();
            nodesG2vars[varShift+keyval.second]->setValue( 0, -vectorVar[0]);
            nodesG2vars[varShift+keyval.second]->setValue( 1,  vectorVar[1]);
            nodesG2vars[varShift+keyval.second]->setValue( 2,  vectorVar[2]);
        }
        
        for ( const auto &keyval : idxs4BoundaryY2fill ) {
            vectorVar = nodesG2vars[varShift+keyval.first]->getValue();
            nodesG2vars[varShift+keyval.second]->setValue( 0,  vectorVar[0]);
            nodesG2vars[varShift+keyval.second]->setValue( 1, -vectorVar[1]);
            nodesG2vars[varShift+keyval.second]->setValue( 2,  vectorVar[2]);
        }
        
        for ( const auto &keyval : idxs4BoundaryZ2fill ) {
            vectorVar = nodesG2vars[varShift+keyval.first]->getValue();
            nodesG2vars[varShift+keyval.second]->setValue( 0,  vectorVar[0]);
            nodesG2vars[varShift+keyval.second]->setValue( 1,  vectorVar[1]);
            nodesG2vars[varShift+keyval.second]->setValue( 2, -vectorVar[2]);
        }
        
        auto end_time = high_resolution_clock::now();
        string msg ="[GridManager] apply BC for var = "+to_string(varName) +"  duration = "
        +to_string(duration_cast<milliseconds>(end_time - start_time).count())+" ms";
        logger->writeMsg(msg.c_str(), DEBUG);
        
        return;
    }
    
    set<int> idealRevertVariables = {ELECTRIC, ELECTRIC_AUX};
    
    if ( idealRevertVariables.count(varName) ){
        for ( const auto &keyval : idxs4BoundaryX2fill ) {
            vectorVar = nodesG2vars[varShift+keyval.first]->getValue();
            nodesG2vars[varShift+keyval.second]->setValue( 0,  vectorVar[0]);
            nodesG2vars[varShift+keyval.second]->setValue( 1, -vectorVar[1]);
            nodesG2vars[varShift+keyval.second]->setValue( 2, -vectorVar[2]);
        }
        for ( const auto &keyval : idxs4BoundaryY2fill ) {
            vectorVar = nodesG2vars[varShift+keyval.first]->getValue();
            nodesG2vars[varShift+keyval.second]->setValue( 0, -vectorVar[0]);
            nodesG2vars[varShift+keyval.second]->setValue( 1,  vectorVar[1]);
            nodesG2vars[varShift+keyval.second]->setValue( 2, -vectorVar[2]);
        }
        for ( const auto &keyval : idxs4BoundaryZ2fill ) {
            vectorVar = nodesG2vars[varShift+keyval.first]->getValue();
            nodesG2vars[varShift+keyval.second]->setValue( 0, -vectorVar[0]);
            nodesG2vars[varShift+keyval.second]->setValue( 1, -vectorVar[1]);
            nodesG2vars[varShift+keyval.second]->setValue( 2,  vectorVar[2]);
        }
        
        auto end_time = high_resolution_clock::now();
        string msg ="[GridManager] apply BC for var = "+to_string(varName) +"  duration = "
        +to_string(duration_cast<milliseconds>(end_time - start_time).count())+" ms";
        logger->writeMsg(msg.c_str(), DEBUG);
        
        return;
    }
    
    // by default set the same value
    for ( const auto &keyval : allBoundaryCells ) {
        vectorVar = nodesG2vars[varShift+keyval.first]->getValue();
        for ( dim = 0; dim < varDim; dim++) {
            nodesG2vars[varShift+keyval.second]->setValue( dim, vectorVar[dim]);
        }
    }
    
    
    auto end_time = high_resolution_clock::now();
    string msg ="[GridManager] apply BC for var = "+to_string(varName) +"  duration = "
    +to_string(duration_cast<milliseconds>(end_time - start_time).count())+" ms";
    logger->writeMsg(msg.c_str(), DEBUG);
    
}


void GridManager::sendRecvIndecis4MPI(){
    int xRes = loader->resolution[0],
        yRes = loader->resolution[1],
        zRes = loader->resolution[2];
    int xResG2 = xRes+2, yResG2 = yRes+2, zResG2 = zRes+2;
    int idx, i, j, k, a, b, c, t;
    vector<int> idxX, idxY, idxZ;
    
    for (a = -1; a <= 1; a++) {
        for (b = -1; b <= 1; b++) {
            for (c = -1; c <= 1; c++) {
                t = (1 + c) + 3 * ((1 + b) + 3 * (1 + a));

                switch (a) {
                    case -1: idxX = {1   , 1   , xRes + 1, xRes + 1}; break;
                    case  0: idxX = {1   , xRes, 1       , xRes    }; break;
                    case  1: idxX = {xRes, xRes, 0       , 0       }; break;
                }
                switch (b) {
                    case -1: idxY = {1   , 1   , yRes + 1, yRes + 1}; break;
                    case  0: idxY = {1   , yRes, 1       , yRes    }; break;
                    case  1: idxY = {yRes, yRes, 0       , 0       }; break;
                }
                switch (c) {
                    case -1: idxZ = {1   , 1   , zRes + 1, zRes + 1}; break;
                    case  0: idxZ = {1   , zRes, 1       , zRes    }; break;
                    case  1: idxZ = {zRes, zRes, 0       , 0       }; break;
                }
                
                counter[t] = (idxX[1] - idxX[0] + 1)*
                             (idxY[1] - idxY[0] + 1)*
                             (idxZ[1] - idxZ[0] + 1);
                sendIdx[t] = new int[counter[t]];
                recvIdx[t] = new int[counter[t]];
                
                int lineIDX = 0;
                for (i = idxX[0]; i <= idxX[1]; i++) {
                    for (j = idxY[0]; j <= idxY[1]; j++) {
                        for (k = idxZ[0]; k <= idxZ[1]; k++) {
                            idx = IDX(i,j,k, xResG2, yResG2, zResG2);
                            sendIdx[t][lineIDX] = idx;
                            lineIDX++;
                        }
                    }
                }

                lineIDX = 0;
                for (i = idxX[2]; i <= idxX[3]; i++) {
                    for (j = idxY[2]; j <= idxY[3]; j++) {
                        for (k = idxZ[2]; k <= idxZ[3]; k++) {
                            idx = IDX(i,j,k, xResG2, yResG2, zResG2);
                            recvIdx[t][lineIDX] = idx;
                            lineIDX++;
                        }
                    }
                }
            }
        }
    }
}


void GridManager::sendBoundary2Neighbor(int varName){
    
    auto start_time = high_resolution_clock::now();
    int varShift = G2nodesNumber*varName;
    int varDim = nodesG2vars[varShift]->getSize();
    int t, i;
    double *sendBuf[27];
    double *recvBuf[27];
    const double* vectorVar;
    
    for (t = 0; t < 27; t++) {
        sendBuf[t] = new double[counter[t]*varDim*sizeof(double)];
        recvBuf[t] = new double[counter[t]*varDim*sizeof(double)];
        for(i = 0; i < counter[t]; i++){
            vectorVar = nodesG2vars[varShift+sendIdx[t][i]]->getValue();
            for (int dim=0;dim<varDim;dim++) {
                sendBuf[t][varDim*i+dim] = vectorVar[dim];
            }
        }
    }
    
    MPI_Status st;
    for (t = 0; t < 27; t++) {
        if (t != 13){
            MPI_Sendrecv(sendBuf[t], counter[t]*varDim, MPI_DOUBLE, loader->neighbors2Send[t], t,
                         recvBuf[t], counter[t]*varDim, MPI_DOUBLE, loader->neighbors2Recv[t], t,
                         MPI_COMM_WORLD, &st);
        }
    }
    
    for (t = 0; t < 27; t++) {
        if (t != 13 && loader->neighbors2Recv[t] != MPI_PROC_NULL) {
            for(i = 0; i < counter[t]; i++){
                for (int dim=0;dim<varDim;dim++) {
                    nodesG2vars[varShift+recvIdx[t][i]]
                    ->setValue( dim, recvBuf[t][varDim*i+dim]);
                }
            }
        }
    }
    
    
    for (t = 0; t < 27; t++) {
        delete [] sendBuf[t];
        delete [] recvBuf[t];
    }
    
    auto end_time = high_resolution_clock::now();
    string msg ="[GridManager] send boundary for var = "+to_string(varName)
                +" duration = "
                +to_string(duration_cast<milliseconds>(end_time - start_time).count())+" ms";
    logger->writeMsg(msg.c_str(), DEBUG);
}



void GridManager::sendRecvIndecis4MPIext(){
    int xRes = loader->resolution[0],
        yRes = loader->resolution[1],
        zRes = loader->resolution[2];
    int xResG2 = xRes+2, yResG2 = yRes+2, zResG2 = zRes+2;
    int idx, i, j, k, a, b, c, t;
    vector<int> idxX, idxY, idxZ;
    
    for (a = -1; a <= 1; a++) {
        for (b = -1; b <= 1; b++) {
            for (c = -1; c <= 1; c++) {
                t = (1 + c) + 3 * ((1 + b) + 3 * (1 + a));
                

                switch (a) {
                    case -1: idxX = {0   , 1       , xRes, xRes + 1}; break;
                    case  0: idxX = {0   , xRes + 1, 0   , xRes + 1}; break;
                    case  1: idxX = {xRes, xRes + 1, 0   , 1       }; break;
                }
                switch (b) {
                    case -1: idxY = {0   , 1       , yRes, yRes + 1}; break;
                    case  0: idxY = {0   , yRes + 1, 0   , yRes + 1}; break;
                    case  1: idxY = {yRes, yRes + 1, 0   , 1       }; break;
                }
                switch (c) {
                    case -1: idxZ = {0   , 1       , zRes, zRes + 1}; break;
                    case  0: idxZ = {0   , zRes + 1, 0   , zRes + 1}; break;
                    case  1: idxZ = {zRes, zRes + 1, 0   , 1       }; break;
                }
                
                counter4Gath[t] = (idxX[1] - idxX[0] + 1)*
                                  (idxY[1] - idxY[0] + 1)*
                                  (idxZ[1] - idxZ[0] + 1);
                sendIdx4Gath[t] = new int[counter4Gath[t]];
                recvIdx4Gath[t] = new int[counter4Gath[t]];
                
                int lineIDX = 0;
                for (i = idxX[0]; i <= idxX[1]; i++) {
                    for (j = idxY[0]; j <= idxY[1]; j++) {
                        for (k = idxZ[0]; k <= idxZ[1]; k++) {
                            idx = IDX(i,j,k, xResG2, yResG2, zResG2);
                            sendIdx4Gath[t][lineIDX] = idx;
                            lineIDX++;
                        }
                    }
                }
                
                lineIDX = 0;
                for (i = idxX[2]; i <= idxX[3]; i++) {
                    for (j = idxY[2]; j <= idxY[3]; j++) {
                        for (k = idxZ[2]; k <= idxZ[3]; k++) {
                            idx = IDX(i,j,k, xResG2, yResG2, zResG2);
                            recvIdx4Gath[t][lineIDX] = idx;
                            lineIDX++;
                        }
                    }
                }
            }
        }
    }
}



void GridManager::gatherBoundaryUsingNeighbor(int varName){
    
    auto start_time = high_resolution_clock::now();
    int varShift = G2nodesNumber*varName;
    int varDim = nodesG2vars[varShift]->getSize();
    
    int t;
    int i, j, k;
    int ijk, ijk0, ijk1;
    double *sendBuf[27];
    double *recvBuf[27];
    MPI_Status st;
    
    int xRes = loader->resolution[0],
        yRes = loader->resolution[1],
        zRes = loader->resolution[2];
    
    const double* varVec;
    
    
    int idx;
    for (t = 0; t < 27; t++) {
        sendBuf[t] = new double[counter4Gath[t]*varDim*sizeof(double)];
        recvBuf[t] = new double[counter4Gath[t]*varDim*sizeof(double)];
        for(i = 0; i<counter4Gath[t]; i++){
            idx = sendIdx4Gath[t][i];
            varVec = nodesG2vars[varShift+idx]->getValue();
            for (int dim=0;dim<varDim;dim++) {
                sendBuf[t][varDim*i+dim] = varVec[dim];
            }
        }
    }
    
    auto end_time1 = high_resolution_clock::now();
    string msg1 ="[GridManager] gatherBoundaryUsingNeighbor: pack for "+to_string(varName)
                +" duration = "+to_string(duration_cast<milliseconds>(end_time1 - start_time).count())+" ms";
    logger->writeMsg(msg1.c_str(), DEBUG);
    
    
    for (t = 0; t < 27; t++) {
        if (t != 13){
            MPI_Sendrecv(sendBuf[t], counter4Gath[t]*varDim,  MPI_DOUBLE, loader->neighbors2Send[t], t,
                         recvBuf[t], counter4Gath[t]*varDim,  MPI_DOUBLE, loader->neighbors2Recv[t], t,
                         MPI_COMM_WORLD, &st);
            
        }
    }
    
    auto end_time2 = high_resolution_clock::now();
    string msg2 ="[GridManager] gatherBoundaryUsingNeighbor: send data for "+to_string(varName)
                +" duration = "+to_string(duration_cast<milliseconds>(end_time2 - end_time1).count())+" ms";
    logger->writeMsg(msg2.c_str(), DEBUG);
    
    for (t = 0; t < 27; t++) {
        if (t != 13 && loader->neighbors2Recv[t] != MPI_PROC_NULL) {
            for(i = 0; i< counter4Gath[t]; i++){
                idx = recvIdx4Gath[t][i];
                for (int dim=0;dim<varDim;dim++) {
                    nodesG2vars[varShift+idx]->addValue( dim, recvBuf[t][varDim*i+dim]);
                }
            }
        }
    }
    
    auto end_time3 = high_resolution_clock::now();
    string msg4 ="[GridManager] gatherBoundaryUsingNeighbor: unpack for "+to_string(varName)
                +" duration = "+to_string(duration_cast<milliseconds>(end_time3 - end_time2).count())+" ms";
    logger->writeMsg(msg4.c_str(), DEBUG);
    
    
    for (t = 0; t < 27; t++) {
        delete [] sendBuf[t];
        delete [] recvBuf[t];
    }
    
    
    if (xRes == 1){
        for (j = 0; j < yRes+2; j++){
            for (k = 0; k < zRes+2; k++){
                ijk0 = IDX(0, j, k, xRes+2, yRes+2, zRes+2);
                ijk  = IDX(1, j, k, xRes+2, yRes+2, zRes+2);
                ijk1 = IDX(2, j, k, xRes+2, yRes+2, zRes+2);
                
                nodesG2vars[varShift+ijk0]
                ->setValue(nodesG2vars[varShift+ijk]->getValue());
                
                nodesG2vars[varShift+ijk1]
                ->setValue(nodesG2vars[varShift+ijk]->getValue());
            }
        }
    }
    
    if (yRes == 1){
        for (i = 0; i < xRes+2; i++){
            for (k = 0; k < zRes+2; k++){
                ijk0 = IDX(i, 0, k, xRes+2, yRes+2, zRes+2);
                ijk  = IDX(i, 1, k, xRes+2, yRes+2, zRes+2);
                ijk1 = IDX(i, 2, k, xRes+2, yRes+2, zRes+2);
                
                nodesG2vars[varShift+ijk0]
                ->setValue(nodesG2vars[varShift+ijk]->getValue());
                
                nodesG2vars[varShift+ijk1]
                ->setValue(nodesG2vars[varShift+ijk]->getValue());
            }
        }
    }
    
    
    if (zRes == 1){
        for (i = 0; i < xRes+2; i++){
            for (j = 0; j < yRes+2; j++){
                ijk0 = IDX(i, j, 0, xRes+2, yRes+2, zRes+2);
                ijk  = IDX(i, j, 1, xRes+2, yRes+2, zRes+2);
                ijk1 = IDX(i, j, 2, xRes+2, yRes+2, zRes+2);
                
                nodesG2vars[varShift+ijk0]
                ->setValue(nodesG2vars[varShift+ijk]->getValue());
                
                nodesG2vars[varShift+ijk1]
                ->setValue(nodesG2vars[varShift+ijk]->getValue());
            }
        }
    }
    
    auto end_time = high_resolution_clock::now();
    string msg ="[GridManager] gatherBoundaryUsingNeighbor: total for "+to_string(varName)
                +" duration = "+to_string(duration_cast<milliseconds>(end_time - start_time).count())+" ms";
    logger->writeMsg(msg.c_str(), DEBUG);

    
}



void GridManager::sendRecvIndecis4MPIonG4(){
    int xRes = loader->resolution[0],
        yRes = loader->resolution[1],
        zRes = loader->resolution[2];
    int xResG2 = xRes+2, yResG2 = yRes+2, zResG2 = zRes+2;
    int xResG4 = xRes+4, yResG4 = yRes+4, zResG4 = zRes+4;
    int idx, i, j, k, a, b, c, t;
    vector<int> idxX, idxY, idxZ;
    
    for (a = -1; a <= 1; a++) {
        for (b = -1; b <= 1; b++) {
            for (c = -1; c <= 1; c++) {
                
                switch (a) {
                    case -1: idxX = {3   , 3       , xRes + 3, xRes + 3}; break;
                    case  0: idxX = {1   , xRes + 2, 1       , xRes + 2}; break;
                    case  1: idxX = {xRes, xRes    , 0       , 0       }; break;
                }
                switch (b) {
                    case -1: idxY = {3   , 3       , yRes + 3, yRes + 3}; break;
                    case  0: idxY = {1   , yRes + 2, 1       , yRes + 2}; break;
                    case  1: idxY = {yRes, yRes    , 0       , 0       }; break;
                }
                switch (c) {
                    case -1: idxZ = {3   , 3       , zRes + 3, zRes + 3}; break;
                    case  0: idxZ = {1   , zRes + 2, 1       , zRes + 2}; break;
                    case  1: idxZ = {zRes, zRes    , 0       , 0       }; break;
                }
                
                t = (1 + c) + 3 * ((1 + b) + 3 * (1 + a));
                
                counterOnG4[t] = (idxX[1] - idxX[0] + 1)*
                                 (idxY[1] - idxY[0] + 1)*
                                 (idxZ[1] - idxZ[0] + 1);
                sendIdxOnG4[t] = new int[counterOnG4[t]];
                recvIdxOnG4[t] = new int[counterOnG4[t]];
                
                int lineIDX = 0;
                for (i = idxX[0]; i <= idxX[1]; i++) {
                    for (j = idxY[0]; j <= idxY[1]; j++) {
                        for (k = idxZ[0]; k <= idxZ[1]; k++) {
                            idx = IDX(i,j,k, xResG4, yResG4, zResG4);
                            sendIdxOnG4[t][lineIDX] = idx;
                            lineIDX++;
                        }
                    }
                }
                lineIDX = 0;
                for (i = idxX[2]; i <= idxX[3]; i++) {
                    for (j = idxY[2]; j <= idxY[3]; j++) {
                        for (k = idxZ[2]; k <= idxZ[3]; k++) {
                            idx = IDX(i,j,k, xResG4, yResG4, zResG4);
                            recvIdxOnG4[t][lineIDX] = idx;
                            lineIDX++;
                        }
                    }
                }
            }
        }
    }
}



void GridManager::smooth(int varName){
    
    int i, j, k, idx, idxG4;
    int xRes = loader->resolution[0],
        yRes = loader->resolution[1],
        zRes = loader->resolution[2];
    int xResG4 = xRes+4, yResG4 = yRes+4, zResG4 = zRes+4;
    int xResG2 = xRes+2, yResG2 = yRes+2, zResG2 = zRes+2;
    int totG4 = xResG4*yResG4*zResG4;
    auto start_time = high_resolution_clock::now();
    int varShift = G2nodesNumber*varName;
    int varDim = nodesG2vars[varShift]->getSize();
    int t;
    double *sendBuf[27];
    double *recvBuf[27];
    const double* vectorVar;
    
    double* varValues = new double[totG4*varDim*sizeof(double)];
    
    for ( i=0; i<xResG2; i++){
        for ( j=0; j<yResG2; j++){
            for ( k=0; k<zResG2; k++){
                idx   = IDX(i  ,j  ,k  ,xResG2,yResG2,zResG2);
                idxG4 = IDX(i+1,j+1,k+1,xResG4,yResG4,zResG4);
                vectorVar = nodesG2vars[varShift+idx]->getValue();
                for (int dim=0;dim<varDim;dim++) {
                    varValues[varDim*idxG4+dim] = vectorVar[dim];
                }
            }
        }
    }
    
    for (t = 0; t < 27; t++) {
        sendBuf[t] = new double[counterOnG4[t]*varDim*sizeof(double)];
        recvBuf[t] = new double[counterOnG4[t]*varDim*sizeof(double)];
        for(i = 0; i < counterOnG4[t]; i++){
            int si = sendIdxOnG4[t][i];
            for (int dim=0;dim<varDim;dim++) {
                sendBuf[t][varDim*i+dim] = varValues[varDim*si+dim];
                
            }
        }
    }
    
    MPI_Status st;
    for (t = 0; t < 27; t++) {
        if (t != 13){
            MPI_Sendrecv(sendBuf[t], counterOnG4[t]*varDim, MPI_DOUBLE, loader->neighbors2Send[t], t,
                         recvBuf[t], counterOnG4[t]*varDim, MPI_DOUBLE, loader->neighbors2Recv[t], t,
                         MPI_COMM_WORLD, &st);
        }
    }
    for (t = 0; t < 27; t++) {
        if (t != 13 && loader->neighbors2Recv[t] != MPI_PROC_NULL) {
            for(i = 0; i < counterOnG4[t]; i++){
                int ri = recvIdxOnG4[t][i];
                for (int dim=0;dim<varDim;dim++) {
                    varValues[varDim*ri+dim] = recvBuf[t][varDim*i+dim];
                }
            }
        }
    }
    
    
    for (t = 0; t < 27; t++) {
        delete [] sendBuf[t];
        delete [] recvBuf[t];
    }
    
    int zeroOrderNeighb[6][3]  =
       {{-1,0 ,0 }, {+1,0 ,0 },
        {0 ,-1,0 }, {0 ,+1,0 },
        {0 ,0 ,-1}, {0 ,0 ,+1}};
    
    int firstOrderNeighb[12][3] =
       {{-1,-1, 0}, {-1,+1, 0},
        {-1, 0,-1}, {-1, 0,+1},
        {0 ,-1,-1}, { 0,-1,+1},
        {0 ,+1,-1}, { 0,+1,+1},
        {+1, 0,-1}, {+1, 0,+1},
        {+1,-1, 0}, {+1,+1, 0}};
    
    int secndOrderNeighb[8][3] =
       {{-1,-1,-1}, {-1,-1,+1},
        {-1,+1,-1}, {-1,+1,+1},
        {+1,-1,-1}, {+1,-1,+1},
        {+1,+1,-1}, {+1,+1,+1}};
    
    const double k2 = 0.125;// 1/8
    const double k3 = 0.0625;// 1/16
    const double k4 = 0.03125;// 1/32
    const double k5 = 0.015625;// 1/64
//                    +--------
//  0.125*1+0.0625*6+0.03125*12+0.015625*8 = 1
    
    int foi, soi, toi;
    int idxFo, idxSo, idxTo;
    
    for ( i=1; i<xResG2+1; i++){
        for ( j=1; j<yResG2+1; j++){
            for ( k=1; k<zResG2+1; k++){
                idxG4 = IDX(i,j,k,xResG4,yResG4,zResG4);
                

                for (int dim=0;dim<varDim;dim++) {
                    
//                    if(abs(varValues[varDim*idxG4+dim]) < EPS8){
//                        continue;
//                    }
                    
                    double smoothedVal = k2*varValues[varDim*idxG4+dim];
                    
                    for(foi = 0; foi<6; foi++){
                        idxFo = IDX(i+zeroOrderNeighb[foi][0],
                                    j+zeroOrderNeighb[foi][1],
                                    k+zeroOrderNeighb[foi][2],
                                    xResG4,yResG4,zResG4);
                        smoothedVal += k3*varValues[varDim*idxFo+dim];
                    }
                    
                    for(soi = 0; soi<12; soi++){
                        idxSo = IDX(i+firstOrderNeighb[soi][0],
                                    j+firstOrderNeighb[soi][1],
                                    k+firstOrderNeighb[soi][2],
                                    xResG4,yResG4,zResG4);
                        smoothedVal += k4*varValues[varDim*idxSo+dim];
                    }
                    
                    for(toi = 0; toi<8; toi++){
                        idxTo = IDX(i+secndOrderNeighb[toi][0],
                                    j+secndOrderNeighb[toi][1],
                                    k+secndOrderNeighb[toi][2],
                                    xResG4,yResG4,zResG4);
                        smoothedVal += k5*varValues[varDim*idxTo+dim];
                    }
                    idx = IDX(i-1,j-1,k-1,xResG2,yResG2,zResG2);
                    nodesG2vars[varShift+idx]->setValue(dim, smoothedVal);
                }
            }
        }
    }
    
    
    delete varValues;
    auto end_time = high_resolution_clock::now();
    string msg ="[GridManager] smooth "+to_string(varName)
                +" duration = "+to_string(duration_cast<milliseconds>(end_time - start_time).count())+" ms";
    logger->writeMsg(msg.c_str(), DEBUG);
    
}

void GridManager::smoothDensAndIonVel(){
    
    int i, j, k, idx, idxG4;
    int xRes = loader->resolution[0],
        yRes = loader->resolution[1],
        zRes = loader->resolution[2];
    int xResG4 = xRes+4, yResG4 = yRes+4, zResG4 = zRes+4;
    int xResG2 = xRes+2, yResG2 = yRes+2, zResG2 = zRes+2;
    int totG4 = xResG4*yResG4*zResG4;
    auto start_time = high_resolution_clock::now();
    int varDim = 4;
    int t;
    double *sendBuf[27];
    double *recvBuf[27];
    const double* vectorVar;
    
    double* densVel = new double[totG4*varDim*sizeof(double)];
    
    for ( i=0; i<xResG2; i++){
        for ( j=0; j<yResG2; j++){
            for ( k=0; k<zResG2; k++){
                idx   = IDX(i  ,j  ,k  ,xResG2,yResG2,zResG2);
                idxG4 = IDX(i+1,j+1,k+1,xResG4,yResG4,zResG4);
                
                vectorVar = nodesG2vars[G2nodesNumber*DENSELEC+idx]->getValue();
                densVel[varDim*idxG4+0] = vectorVar[0];
                vectorVar = nodesG2vars[G2nodesNumber*VELOCION+idx]->getValue();
                for (int dim=1;dim<varDim;dim++) {
                    densVel[varDim*idxG4+dim] = vectorVar[dim-1];
                }
            }
        }
    }
    
    for (t = 0; t < 27; t++) {
        sendBuf[t] = new double[counterOnG4[t]*varDim*sizeof(double)];
        recvBuf[t] = new double[counterOnG4[t]*varDim*sizeof(double)];
        for(i = 0; i < counterOnG4[t]; i++){
            int si = sendIdxOnG4[t][i];
            for (int dim=0;dim<varDim;dim++) {
                sendBuf[t][varDim*i+dim] = densVel[varDim*si+dim];
                
            }
        }
    }
    
    MPI_Status st;
    for (t = 0; t < 27; t++) {
        if (t != 13){
            MPI_Sendrecv(sendBuf[t], counterOnG4[t]*varDim, MPI_DOUBLE, loader->neighbors2Send[t], t,
                         recvBuf[t], counterOnG4[t]*varDim, MPI_DOUBLE, loader->neighbors2Recv[t], t,
                         MPI_COMM_WORLD, &st);
        }
    }
    for (t = 0; t < 27; t++) {
        if (t != 13 && loader->neighbors2Recv[t] != MPI_PROC_NULL) {
            for(i = 0; i < counterOnG4[t]; i++){
                int ri = recvIdxOnG4[t][i];
                for (int dim=0;dim<varDim;dim++) {
                    densVel[varDim*ri+dim] = recvBuf[t][varDim*i+dim];
                }
            }
        }
    }
    
    
    for (t = 0; t < 27; t++) {
        delete [] sendBuf[t];
        delete [] recvBuf[t];
    }
    
    
    
    
    
    /*
     *                           |               |               |
     *                                           |
     *                  (i-1,j+1).B______________A.(i,j+1)________.B(i+1,j+1)
     *                           .               |               .
     *                           |               |               |
     *                           .               |               .
     *                                           |               |
     *                           .               |               .
     *                           |               |               |
     *                           .               |               .
     *                           |               |               |
     *                   A(i-1,j) ._______________O.(i,j)___________.A(i+1,j)
     *                           |               |               |
     *                           .               |               .
     *                           |               |               |
     *                                           |               .
     *                           |               |               |
     *                           .               |               .
     *                           |               |               |
     *                           .               |               .
     *                  (i-1,j-1)B.______________A.(i,j-1)________.B(i+1,j-1)
     *                               
     */


    int zeroOrderNeighb[6][3]  =
       {{-1,0 ,0 }, {+1,0 ,0 },
        {0 ,-1,0 }, {0 ,+1,0 },
        {0 ,0 ,-1}, {0 ,0 ,+1}};
    
    int firstOrderNeighb[12][3] =
       {{-1,-1, 0}, {-1,+1, 0},
        {-1,0 ,-1}, {-1,0 ,+1},
        {0 ,-1,-1}, {0 ,-1,+1},
        {0 ,+1,-1}, {0 ,+1,+1},
        {+1,0 ,-1}, {+1,0 ,+1},
        {+1,-1, 0}, {+1,+1, 0}};
    
    int secndOrderNeighb[8][3] =
       {{-1,-1,-1}, {-1,-1,+1},
        {-1,+1,-1}, {-1,+1,+1},
        {+1,-1,-1}, {+1,-1,+1},
        {+1,+1,-1}, {+1,+1,+1}};
    
    const double k2 = 0.125;
    const double k3 = 0.0625;
    const double k4 = 0.03125;
    const double k5 = 0.015625;
    
    int foi, soi, toi;
    int idxFo, idxSo, idxTo;
    
    for ( i=1; i<xResG2+1; i++){
        for ( j=1; j<yResG2+1; j++){
            for ( k=1; k<zResG2+1; k++){
                idxG4 = IDX(i,j,k,xResG4,yResG4,zResG4);
                
                double smoothedVal[4] = {0.0,0.0,0.0,0.0};
                
                for (int dim=0;dim<varDim;dim++) {
                    
                    if(abs(densVel[varDim*idxG4+dim]) < EPS8){
                        continue;
                    }
                    
                    smoothedVal[dim] = k2*densVel[varDim*idxG4+dim];
                    
                    for(foi = 0; foi<6; foi++){
                        idxFo = IDX(i+zeroOrderNeighb[foi][0],
                                    j+zeroOrderNeighb[foi][1],
                                    k+zeroOrderNeighb[foi][2],
                                    xResG4,yResG4,zResG4);
                        smoothedVal[dim] += k3*densVel[varDim*idxFo+dim];
                    }
                
                    for(soi = 0; soi<12; soi++){
                        idxSo = IDX(i+firstOrderNeighb[soi][0],
                                    j+firstOrderNeighb[soi][1],
                                    k+firstOrderNeighb[soi][2],
                                    xResG4,yResG4,zResG4);
                        smoothedVal[dim] += k4*densVel[varDim*idxSo+dim];
                    }
               
                    for(toi = 0; toi<8; toi++){
                        idxTo = IDX(i+secndOrderNeighb[toi][0],
                                    j+secndOrderNeighb[toi][1],
                                    k+secndOrderNeighb[toi][2],
                                    xResG4,yResG4,zResG4);
                        smoothedVal[dim] += k5*densVel[varDim*idxTo+dim];
                    }
                    
                }
                idx = IDX(i-1,j-1,k-1,xResG2,yResG2,zResG2);
                nodesG2vars[G2nodesNumber*DENSELEC+idx]->setValue(0, smoothedVal[0]);
                for (int dim=1;dim<varDim;dim++) {
                    nodesG2vars[G2nodesNumber*VELOCION+idx]->setValue(dim-1, smoothedVal[dim]);
                }
            }
        }
    }
    
    
    delete [] densVel;
    auto end_time = high_resolution_clock::now();
    string msg ="[GridManager] smooth N V duration = "
                +to_string(duration_cast<milliseconds>(end_time - start_time).count())+" ms";
    logger->writeMsg(msg.c_str(), DEBUG);

}

vector<vector<VectorVar>> GridManager::getVectorVariablesForAllNodes(){
    
    int idxG2, idxG1;
    int xRes = loader->resolution[0], yRes = loader->resolution[1], zRes = loader->resolution[2];
    int xResG1 = xRes+1, yResG1 = yRes+1, zResG1 = zRes+1;
    int xResG2 = xRes+2, yResG2 = yRes+2, zResG2 = zRes+2;
    vector<vector<VectorVar>> result;
    result.reserve(xRes*yRes*zRes);
    int i,j,k;
    
    set<int> stopList = {ELECTRIC_AUX, CURRENT, VELOCELE, DRIVER, DRIVER_AUX, PRESSURE_AUX, PRESSURE_SMO};
    int numOfSpecies = loader->getNumberOfSpecies();
    for(i=0;i<numOfSpecies;i++){
        stopList.insert(DENS_AUX(i));
    }
    for ( i=0; i<xRes; i++){
        for ( j=0; j<yRes; j++){
            for ( k=0; k<zRes; k++){
                idxG1 = IDX(i  ,j  ,k  ,xResG1,yResG1,zResG1);
                idxG2 = IDX(i+1,j+1,k+1,xResG2,yResG2,zResG2);
                vector<VectorVar> allVars;
                for ( int varN=0; varN<totVarsOnG2; varN++){
                    if(stopList.count(varN)){
                        continue;
                    }
                    allVars.push_back(*nodesG2vars[G2nodesNumber*varN+idxG2]);
                }
                allVars.push_back(*nodesG1vars[G1nodesNumber*MAGNETIC+idxG1]);
                result.push_back(allVars);
            }
        }
    }
    
    return result;
}



void GridManager::setVectorVariableForNodeG1(int idx, VectorVar variable){
    nodesG1vars[G1nodesNumber*variable.getName()+idx]->setValue(variable.getValue());
}


void GridManager::setVectorVariableForNodeG2(int idx, VectorVar variable){
    nodesG2vars[G2nodesNumber*variable.getName()+idx]->setValue(variable.getValue());
}

void GridManager::setVectorVariableForNodeG2(int idx, int name, int dim, double value){
    nodesG2vars[G2nodesNumber*name+idx]->setValue(dim, value);
}

void GridManager::setVectorVariableForNodeG1(int idx, int name, int dim, double value){
    nodesG1vars[G1nodesNumber*name+idx]->setValue(dim, value);
}

void GridManager::addVectorVariableForNodeG2(int idx, int name, int dim, double value){
    nodesG2vars[G2nodesNumber*name+idx]->addValue(dim, value);
}


VectorVar** GridManager::getVectorVariableOnG1(int varName){
    return &nodesG1vars[G1nodesNumber*varName];
}


VectorVar** GridManager::getVectorVariableOnG2(int varName){
     return &nodesG2vars[G2nodesNumber*varName];
}

