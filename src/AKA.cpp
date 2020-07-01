#include "core/SimulationManager.hpp"

#include <mpi.h>

int main(int ac, char **av) {

    //init MPI
    MPI_Init(&ac, &av);
    
    SimulationManager simMng(ac, av);
    simMng.runSimulation(ac, av);

    //close MPI
    sleep(1);
    MPI_Finalize();
    return 0;
}
