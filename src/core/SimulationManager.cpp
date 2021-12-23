#include "SimulationManager.hpp"
#include <iostream>
#include <string>
#include <thread>
#include <unistd.h>
using namespace std;
using namespace chrono;


SimulationManager::SimulationManager(int ac, char **av) {
    
    logger.reset(new Logger());
    
    int rank ;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    //init random
    srand(time(NULL) + rank);
    
    unsigned int pause = 20000;
    usleep(rank*pause);

    
    if(rank == 0){
        char buff[20];
        time_t now = time(NULL);
        strftime(buff, 20, "%Y-%m-%d %H:%M:%S", localtime(&now));
        unsigned int fancyStart = 0;//300000;
        logger->writeMsg("***********************************************************", INFO);
        usleep(fancyStart);
        logger->writeMsg("****                                                   ****", INFO);
        usleep(fancyStart);
        logger->writeMsg(string("****     START  "+string(buff)+"                    ****").c_str(), INFO);
        usleep(fancyStart);
        logger->writeMsg("****                                                   ****", INFO);
        usleep(fancyStart);
        logger->writeMsg("****       O        0  A       A         55555  2222   ****", INFO);
        usleep(fancyStart);
        logger->writeMsg("****      O O       O  A      A A        5          2  ****", INFO);
        usleep(fancyStart);
        logger->writeMsg("****     O   O      O A      A   A       5          2  ****", INFO);
        usleep(fancyStart);
        logger->writeMsg("****    OOOOOOO     0A      AAAAAAA      5555    2222  ****", INFO);
        usleep(fancyStart);
        logger->writeMsg("****   O       0    O A    A       A         5  2      ****", INFO);
        usleep(fancyStart);
        logger->writeMsg("****  O         O   O  A  A         A    55555  22222  ****", INFO);
        usleep(fancyStart);
        logger->writeMsg("****   Arbitrary  Kinetic  Algorithm                   ****", INFO);
        usleep(fancyStart);
        logger->writeMsg("***********************************************************", INFO);
        usleep(fancyStart);
        
    }
    
    
    
    if (ac > 1){

        setenv("INPUTFILEPATH", av[1] , 1);
	auto pythonPath = getenv("PYTHONPATH");
	string msg;
	string pyPath;
	if (pythonPath == NULL){
		msg = "[SimulationManager] variable PYTHONPATH is not set !!!";
		pyPath = av[1];
	}else{
		msg = "[SimulationManager] variable PYTHONPATH ="+string(pythonPath);
		pyPath = string(pythonPath) +":"+av[1];
	}	 

	logger->writeMsg(msg.c_str(), INFO);
	setenv("PYTHONPATH", pyPath.c_str() , 1);
    }else{
        if(rank == 0){
            string msg ="[SimulationManager] ATTENTION!!! used default input file path src/input/ !!!";
            logger->writeMsg(msg.c_str(), CRITICAL);
        }
        const char  *PRJ_PATH = "src/input/";
        setenv("INPUTFILEPATH", PRJ_PATH , 1);
    }
    initialize();
}


void SimulationManager::initialize() {
    

    loader.reset(new Loader());
    loader->load();
        
    gridMng.reset(new GridManager(loader));
    boundMng.reset(new BoundaryManager(loader));
    
    pusher.reset(new Pusher(loader, gridMng, boundMng));
    
    initMng.reset(new ModelInitializer(loader, gridMng, pusher));
    
    laserMng.reset(new LaserMockManager(loader, gridMng, pusher));
    hydroMng.reset(new HydroManager(loader, gridMng, pusher));
    closureMng.reset(new ClosureManager(loader, gridMng));
    elemagMng.reset(new EleMagManager(loader, gridMng));
    
    solver.reset(new Solver(loader, gridMng, pusher,
                                  hydroMng, elemagMng, closureMng,
                                  laserMng));    
    
    writer.reset(new Writer(loader, gridMng, pusher));
    
    logger->writeMsg("[SimulationManager] init...OK", DEBUG);
}

void SimulationManager::runSimulation(int ac, char **av) {
    logger->writeMsg("[SimulationManager] launch simulation...OK", DEBUG);
    auto start_time_tot = high_resolution_clock::now();
    
    int maxTimeStep = loader->getMaxTimestepsNum();
    int maxTimeStep2Write = loader->getTimestepsNum2Write();
    int fileNumCount = 0;
    int i_time;
    int STOP_SIMULATION = 1;

    int rank ;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  
    for ( i_time=0; i_time < maxTimeStep; i_time++ ){
        auto start_time = high_resolution_clock::now();
        if(rank == 0){
            logger->writeMsg("*****************************************************", INFO);
            logger->writeMsg("****                                             ****", INFO);
            logger->writeMsg(string("[SimulationManager] step = "
                                    +to_string(i_time)).c_str(), INFO);
            logger->writeMsg(string("[SimulationManager] time = "
                                    +to_string(i_time*loader->getTimeStep())).c_str(), INFO);
            logger->writeMsg("****                                             ****", INFO);
            logger->writeMsg("*****************************************************", INFO);
        }
        if( i_time % maxTimeStep2Write == 0 ){
            #ifdef GET_ION_PRESSURE
            hydroMng->setIonPressureTensor();
            #endif

            writer->write(fileNumCount);
//            writer->writeAllForRestart();
            fileNumCount++;
        }
        

        if(STOP_SIMULATION == 0){
             break;
        }
        
        
        if( solver->solve(i_time) == SOLVE_FAIL ){
            STOP_SIMULATION = 0;
            
            if(rank == 0){
                logger->writeMsg("*****************************************************", CRITICAL);
                logger->writeMsg("****                                             ****", CRITICAL);
                logger->writeMsg("****  STOP SIMULATION!!!    UNEXPECTED PROBLEM   ****", CRITICAL);
                logger->writeMsg("****                                             ****", CRITICAL);
                logger->writeMsg("****   try debug mode 'make -DLOG'               ****", CRITICAL);
                logger->writeMsg("*****************************************************", CRITICAL);
                exit(-1);
            }
        }
        
        if(rank == 0){
            auto end_time = high_resolution_clock::now();
            string msg = "[SimulationManager] Step duration = "
                        +to_string(duration_cast<milliseconds>(end_time - start_time).count())+" ms";
            logger->writeMsg(msg.c_str(), INFO);
            logger->writeMsg("xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx", INFO);
        }
    }
    
    if(rank == 0){
        auto end_time_tot = high_resolution_clock::now();
        logger->writeMsg("*****************************************************", INFO);
        logger->writeMsg("****                                             ****", INFO);
        logger->writeMsg("****          FINISH SIMULATION                  ****", INFO);
        string msg = "**** Total duration for "+to_string(i_time)+" steps : "
                    +to_string(duration_cast<minutes>(end_time_tot - start_time_tot).count())+" min";
        logger->writeMsg(msg.c_str(), INFO);
        logger->writeMsg("****                                             ****", INFO);
        logger->writeMsg("*****************************************************", INFO);
    }

}



void SimulationManager::finilize() {
    logger->writeMsg("[SimulationManager] finalize...OK", DEBUG);
}

SimulationManager::~SimulationManager() {
    finilize();
}

