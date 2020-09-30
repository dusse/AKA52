#include "Solver.hpp"



using namespace std;
using namespace chrono;




Solver::Solver(shared_ptr<Loader> load,
                           shared_ptr<GridManager> grid,
                           shared_ptr<Pusher> pshr,
                           shared_ptr<HydroManager> hydro,
                           shared_ptr<EleMagManager> em,
                           shared_ptr<ClosureManager> closure,
                           shared_ptr<LaserMockManager> lm):loader(move(load)),
                            gridMng(move(grid)), pusher(move(pshr)),
                            hydroMng(move(hydro)), emMng(move(em)),
                            closureMng(move(closure)), laserMng(move(lm))
{
    logger.reset(new Logger());
    
    initialize();
    logger->writeMsg("[Solver] create...OK", DEBUG);
}

void Solver::initialize()
{
    logger->writeMsg("[Solver] initialize...OK", DEBUG);
}


int Solver::solve(int i_time)
{
    int rank ;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    
    try{
        performCalculation(PREDICTOR, i_time);
        performCalculation(CORRECTOR, i_time);
        
    }catch(...){
        return SOLVE_FAIL;
    }

    if (rank == 0){
        double ts = loader->getTimeStep();
        string msg ="[Solver] SOLVER step ="
                        +to_string(i_time)+"; time = "+to_string(i_time*ts);
        logger->writeMsg(msg.c_str(), DEBUG);
    }
    
    return SOLVE_OK;
}

void Solver::performCalculation(int PHASE, int i_time){
    logger->writeMsg("[Solver] performCaclculation...", DEBUG);
    
    pusher->push(PHASE, i_time);
    
    if ( loader->numOfSpots > 0 ) {
        laserMng->addIons();
    }
    
    hydroMng->gatherMoments(PHASE);

    emMng->calculateBhalf(PHASE);
    emMng->calculateJhalf(PHASE);
    
    closureMng->calculatePressure(PHASE, i_time);
    
    if ( loader->numOfSpots > 0 ) {
        laserMng->accelerate(i_time);
    }
    
    emMng->calculateEnext(PHASE);

    emMng->calculateBnext();
    emMng->calculateJnext();
    
    logger->writeMsg("[Solver] performCaclculation...OK", DEBUG);
}



Solver::~Solver(){
    finilize();
}

void Solver::finilize(){
    
}
