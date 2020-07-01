#ifndef GridManager_hpp
#define GridManager_hpp
#include <stdio.h>
#include <vector>
#include <map>
#include <set>
#include <cmath>
#include <memory>
#include <chrono>
#include <mpi.h>
#include "../misc/Logger.hpp"
#include "../misc/Misc.hpp"
#include "../input/Loader.hpp"

#include "../common/variables/VectorVar.hpp"


enum G1VAR{
    MAGNETIC,
    MAGNETIC_AUX,
    SIZEG1
};

enum G2VAR{
    ELECTRIC,
    ELECTRIC_AUX,
    CURRENT,
    CURRENT_AUX,
    VELOCION,
    DENSELEC,
    PRESSURE,
    PRESSURE_AUX,
    DRIVER,
    DRIVER_AUX,
    SIZEG2
};


class GridManager{
    
    
private:
    
    std::unique_ptr<Logger> logger;
    std::shared_ptr<Loader> loader;
    
    
    int NUM_OF_MAIN_G2VARS = SIZEG2;
    int SHIFT_MAIN_DENS = 0;
    int SHIFT_MAIN_DENS_AUX  = 0;
    
    
    int G1nodesNumber;
    int G2nodesNumber;
    int totVarsOnG2;
    int totVarsOnG1;
    int totalNodeNumber;
    
    VectorVar** nodesG1vars;
    VectorVar** nodesG2vars;
    
    int* neibors4G2spatialDerX;
    int* neibors4G2spatialDerY;
    int* neibors4G2spatialDerZ;
    
    int* negbors4PresX;
    int* negbors4PresY;
    int* negbors4PresZ;
    
    int* negbors4LaplacX;
    int* negbors4LaplacY;
    int* negbors4LaplacZ;
    
    int* neibors4G1spatialDerX;
    int* neibors4G1spatialDerY;
    int* neibors4G1spatialDerZ;
    
    std::map<int, int> idxs4BoundaryX2fill;
    std::map<int, int> idxs4BoundaryY2fill;
    std::map<int, int> idxs4BoundaryZ2fill;
    
    int* neighbourhood;
    
    int counter[27], counter4Gath[27], counterOnG4[27];
    
    int* sendIdx[27];
    int* recvIdx[27];
    
    int* sendIdx4Gath[27];
    int* recvIdx4Gath[27];
    
    int* sendIdxOnG4[27];
    int* recvIdxOnG4[27];
    
    void initialize();
    
    void initG1Nodes();
    void initG2Nodes();
    
    void sendRecvIndecis4MPI();
    void sendRecvIndecis4MPIext();
    void sendRecvIndecis4MPIonG4();
    
    void initBoundaryIndecies();
    
public:
    
    GridManager(std::shared_ptr<Loader>);
    ~GridManager();
   
    int DENS_AUX(int);
    int DENS_VEL(int);
    
    int getVarsNumOnG2();
    int getVarsNumOnG1();
    const int* getNghbd4LaplacInXonG2();
    const int* getNghbd4LaplacInYonG2();
    const int* getNghbd4LaplacInZonG2();
    
    const int* getNghbd4DivPeInXonG2();
    const int* getNghbd4DivPeInYonG2();
    const int* getNghbd4DivPeInZonG2();
    
    const int* getNeibors4G2spatialDerX();
    const int* getNeibors4G2spatialDerY();
    const int* getNeibors4G2spatialDerZ();
    
    const int* getNeibors4G1spatialDerX();
    const int* getNeibors4G1spatialDerY();
    const int* getNeibors4G1spatialDerZ();
    
    VectorVar** getVectorVariableOnG1(int);
    VectorVar** getVectorVariableOnG2(int);
    std::vector<std::vector<VectorVar>> getVectorVariablesForAllNodes();
    
    void setVectorVariableForNodeG1(int, VectorVar);
    void setVectorVariableForNodeG2(int, VectorVar);
    
    void setVectorVariableForNode(int, VectorVar);
    
    void setVectorVariableForNodeG1(int , int , int, double);
    void setVectorVariableForNodeG2(int , int , int, double);
    void addVectorVariableForNodeG2(int , int , int, double);
    
    const int* getNeighbourhoodOnG1();
    
    void sendBoundary2Neighbor(int);
    void gatherBoundaryUsingNeighbor(int);
    void applyBC(int);
    
    void smoothDensAndIonVel();
    void smooth(int);
    
};
#endif /* GridManager_hpp */
