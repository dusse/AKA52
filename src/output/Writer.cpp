#include "Writer.hpp"

using namespace std;
using namespace chrono;

Writer::Writer(shared_ptr<Loader> load,
               shared_ptr<GridManager> gridMnr,
               shared_ptr<Pusher> push):
                loader(move(load)),
                gridMgr(move(gridMnr)),
                pusher(move(push)){
                    
    logger.reset(new Logger());

    this->outputDir       = loader->getOutputDir();
    this->fileNamePattern = loader->getFilenameTemplate();
    
    logger->writeMsg("[Writer] initialize...OK", DEBUG);
}

void Writer::write(int fileNum){
    auto start_time = high_resolution_clock::now();
   
    logger->writeMsg("[Writer] writing...", DEBUG);
    int subdomainXSize = loader->resolution[0];
    int subdomainYSize = loader->resolution[1];
    int subdomainZSize = loader->resolution[2];
    
    int size[3] = {subdomainXSize, subdomainYSize, subdomainZSize};
    
    int totalNodeNum = subdomainXSize*subdomainYSize*subdomainZSize;
    int ijNode;
    vector<vector<VectorVar>> vectorVars = gridMgr->getVectorVariablesForAllNodes();
    vector<string> humanReadableNames = gridMgr->getHumanReadableOutputVarNames();
    
    MPI_Info info = MPI_INFO_NULL;
    
    hid_t access = H5Pcreate(H5P_FILE_ACCESS);
    H5Pset_fapl_mpio(access, MPI_COMM_WORLD, info);
    
    string fileName = outputDir + fileNamePattern + to_string(fileNum) + ".h5";
    hid_t fileID = H5Fcreate(fileName.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, access);
    
    const string groupname = "/vars";
    hid_t group   = H5Gcreate(fileID, groupname.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    
    hid_t dxpl_id = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(dxpl_id, H5FD_MPIO_COLLECTIVE);
    
    int vars4EachNode = vectorVars[0].size();
    for (int idx_var = 0; idx_var < vars4EachNode; idx_var++) {
        
        int varDim = vectorVars[0][idx_var].getSize();
        for (int dir = 0; dir < varDim; dir++) {
            
	    string varName  = humanReadableNames[idx_var] +"_"+ to_string(dir);
            double* var = new double[totalNodeNum];
            
            for (ijNode = 0; ijNode < totalNodeNum; ijNode++) {
                var[ijNode] = (vectorVars[ijNode][idx_var].getValue())[dir];
            }
            
            writeParallel(fileID, group, dxpl_id, varName, var, size);
            
            delete [] var;
        }
    }
    
    H5Gclose(group);
    H5Fflush(fileID, H5F_SCOPE_GLOBAL);
    H5Pclose(dxpl_id);
    H5Fclose(fileID);
    
    auto end_time = high_resolution_clock::now();
    string msg ="[Writer] writing duration = "
                +to_string(duration_cast<milliseconds>(end_time - start_time).count())+" ms";
    logger->writeMsg(msg.c_str(), DEBUG);
}




void Writer::writeParallel(hid_t fileID, hid_t group, hid_t dxpl_id,
                           string dsetName , const double* data, int sizes[3]){

    hsize_t offset[3];
    hsize_t nLoc[3];
    hsize_t nGlob[3];

    for (int dir = 0; dir < 3; dir++) {
        offset[dir] = loader->offsetInPixels[dir];
        nLoc[dir]   = sizes[dir];
        nGlob[dir]  = loader->totPixelsPerBoxSide[dir];
    }

    hid_t memspace  = H5Screate_simple(3, nLoc , NULL);
    hid_t filespace = H5Screate_simple(3, nGlob, NULL);
    H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset, NULL, nLoc, NULL);

    hid_t dset;
    dset  = H5Dcreate(group, dsetName.c_str(), H5T_NATIVE_DOUBLE, filespace,
                      H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    
    H5Dwrite(dset, H5T_NATIVE_DOUBLE, memspace, filespace, dxpl_id, data);
    H5Dclose(dset);

    H5Sclose(memspace);
    H5Sclose(filespace);
}



void Writer::writeParticles(int fileNum){
    auto start_time = high_resolution_clock::now();
    logger->writeMsg("[Writer] writing particles ...", DEBUG);
    int rank ;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Info info = MPI_INFO_NULL;
    
    int coresNum;
    MPI_Comm_size(MPI_COMM_WORLD, &coresNum);
    int* particlesPerCore = new int[coresNum];
    Particle** particles = pusher->getParticles();
    int totalPrtclNumber = pusher->getTotalParticleNumber();
    int totalUnfrozenPrtclNumber = 0;
    int idx, type;
    for( idx = 0; idx < totalPrtclNumber; idx++ ){
        type = particles[idx]->getType();
        if( loader->getIfSpeciesFrozen(type) == 0 ){
            totalUnfrozenPrtclNumber++;
        }
    }
    MPI_Allgather(&totalUnfrozenPrtclNumber, 1, MPI_INT, particlesPerCore, 1, MPI_INT, MPI_COMM_WORLD);

    int partoffset = 0, partglob = 0;
    for( int coreNum = 0; coreNum < coresNum; coreNum++ ){
        partglob += particlesPerCore[coreNum];
        if( coreNum < rank ){
            partoffset += particlesPerCore[coreNum];
        }
    }
    delete[] particlesPerCore;
    
    hid_t access = H5Pcreate(H5P_FILE_ACCESS);
    H5Pset_fapl_mpio(access, MPI_COMM_WORLD, info);
    
    string fileName = outputDir + fileNamePattern +  "particles_"+ to_string(fileNum)+".h5";
    hid_t fileID = H5Fcreate(fileName.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, access);
    
    const string groupname = "/vars";
    hid_t group   = H5Gcreate(fileID, groupname.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    
    hid_t dxpl_id = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(dxpl_id, H5FD_MPIO_COLLECTIVE);
 
    const int SHORT_PARTICLE_SIZE = 7;
    double* particles2save = new double[SHORT_PARTICLE_SIZE*totalUnfrozenPrtclNumber];
    double* prtclPos;
    double* prtclVel;
    int ind2Write = 0;
    for( idx = 0; idx < totalPrtclNumber; idx++ ){
            prtclPos = particles[idx]->getPosition();
            prtclVel = particles[idx]->getVelocity();
            type     = particles[idx]->getType();
            if( loader->getIfSpeciesFrozen(type) == 0 ){
                int shift = SHORT_PARTICLE_SIZE*ind2Write;
                particles2save[shift+0] = type;
                particles2save[shift+1] = 0.5*(prtclPos[0]+prtclPos[3]);
                particles2save[shift+2] = 0.5*(prtclPos[1]+prtclPos[4]);
                particles2save[shift+3] = 0.5*(prtclPos[2]+prtclPos[5]);
                particles2save[shift+4] = 0.5*(prtclVel[0]+prtclVel[3]);
                particles2save[shift+5] = 0.5*(prtclVel[1]+prtclVel[4]);
                particles2save[shift+6] = 0.5*(prtclVel[2]+prtclVel[5]);
                ind2Write++;
            }
    }

    writeParallelWithOffset(fileID, group, dxpl_id, "particles", particles2save,
                            SHORT_PARTICLE_SIZE*totalUnfrozenPrtclNumber,
                            SHORT_PARTICLE_SIZE*partglob,
                            SHORT_PARTICLE_SIZE*partoffset);

    
    H5Fflush(fileID, H5F_SCOPE_GLOBAL);
    delete[] particles2save;    
    H5Gclose(group);
    H5Pclose(dxpl_id);   
    H5Fclose(fileID);    
    
    auto end_time = high_resolution_clock::now();
    string msg ="[Writer] writing particles duration = "
    +to_string(duration_cast<milliseconds>(end_time - start_time).count())+" ms";
    logger->writeMsg(msg.c_str(), DEBUG);
}

void Writer::writeLeftParticles(int fileNum){
    auto start_time = high_resolution_clock::now();
    // logger->writeMsg("[Writer] writing left particles ...", DEBUG);
    int idx, type;
    int rank, coresNum ;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Info info = MPI_INFO_NULL;    
    MPI_Comm_size(MPI_COMM_WORLD, &coresNum);

    int* particlesPerCore = new int[coresNum];
    vector<shared_ptr<Particle>> particles = pusher->getLeftParticles();
    int totalPrtclNumber = particles.size();
    logger->writeMsg(("[Writer] writing left particles ... total left particles number = "
        +to_string(totalPrtclNumber)).c_str(), DEBUG);
    int totalLeftPrtclNumber = totalPrtclNumber;
    MPI_Allgather(&totalLeftPrtclNumber, 1, MPI_INT, particlesPerCore, 1, MPI_INT, MPI_COMM_WORLD);

    int partoffset = 0, partglob = 0;
    for( int coreNum = 0; coreNum < coresNum; coreNum++ ){
        partglob += particlesPerCore[coreNum];
        if( coreNum < rank ){
            partoffset += particlesPerCore[coreNum];
        }
    }
    delete[] particlesPerCore;
    
    hid_t access = H5Pcreate(H5P_FILE_ACCESS);
    H5Pset_fapl_mpio(access, MPI_COMM_WORLD, info);
    
    string fileName = outputDir + fileNamePattern +  "left_particles_"+ to_string(fileNum)+".h5";
    hid_t fileID = H5Fcreate(fileName.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, access);
    
    const string groupname = "/vars";
    hid_t group   = H5Gcreate(fileID, groupname.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    
    hid_t dxpl_id = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(dxpl_id, H5FD_MPIO_COLLECTIVE);
 
    const int SHORT_PARTICLE_SIZE = 7;
    double* particles2save = new double[SHORT_PARTICLE_SIZE*totalLeftPrtclNumber];
    double* prtclPos;
    double* prtclVel;
    int ind2Write = 0;
    for( idx = 0; idx < totalPrtclNumber; idx++ ){
            prtclPos = particles[idx]->getPosition();
            prtclVel = particles[idx]->getVelocity();
            type     = particles[idx]->getType();
            int shift = SHORT_PARTICLE_SIZE*idx;
            particles2save[shift+0] = type;
            particles2save[shift+1] = prtclPos[0];
            particles2save[shift+2] = prtclPos[1];
            particles2save[shift+3] = prtclPos[2];
            particles2save[shift+4] = prtclVel[0];
            particles2save[shift+5] = prtclVel[1];
            particles2save[shift+6] = prtclVel[2];            
    }
    
    writeParallelWithOffset(fileID, group, dxpl_id, "particles", particles2save,
                            SHORT_PARTICLE_SIZE*totalLeftPrtclNumber,
                            SHORT_PARTICLE_SIZE*partglob,
                            SHORT_PARTICLE_SIZE*partoffset);

    
    H5Fflush(fileID, H5F_SCOPE_GLOBAL);
    delete[] particles2save;    
    H5Gclose(group);
    H5Pclose(dxpl_id);   
    H5Fclose(fileID);    
    
    auto end_time = high_resolution_clock::now();
    string msg ="[Writer] writing left particles duration = "
    +to_string(duration_cast<milliseconds>(end_time - start_time).count())+" ms";
    logger->writeMsg(msg.c_str(), DEBUG);
}


void Writer::writeParallelWithOffset(hid_t fileID, hid_t group, hid_t dxpl_id,
                                   string dsetName , const double* data,
                                     int sizeLoc, int sizeGlob, int off){
    
    hsize_t offset[1];
    hsize_t nLoc[1];
    hsize_t nGlob[1];
    offset[0] = off;
    nLoc[0] = sizeLoc;
    nGlob[0] = sizeGlob;
    
    hid_t memspace  = H5Screate_simple(1, nLoc , NULL);
    hid_t filespace = H5Screate_simple(1, nGlob, NULL);
    H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset, NULL, nLoc, NULL);
    
    hid_t dset;
    dset  = H5Dcreate(group, dsetName.c_str(), H5T_NATIVE_DOUBLE, filespace,
                      H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    
    H5Dwrite(dset, H5T_NATIVE_DOUBLE, memspace, filespace, dxpl_id, data);
    
    H5Dclose(dset);
    
    H5Sclose(memspace);
    H5Sclose(filespace);
}


void Writer::writeAttributeInt(hid_t fileID, hid_t group, hid_t dxpl_id,
                    string name , int val){
    int rank ;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    
    int coresNum;
    MPI_Comm_size(MPI_COMM_WORLD, &coresNum);
    
    hsize_t offset, nLoc, nGlob;
    offset = rank;
    nLoc = 1;
    nGlob = coresNum;
    
    int * data = new int[1];
    data[0] = val;
    hid_t memspace  = H5Screate_simple(1, &nLoc , NULL);
    hid_t filespace = H5Screate_simple(1, &nGlob, NULL);
    H5Sselect_hyperslab(filespace, H5S_SELECT_SET, &offset, NULL, &nLoc, NULL);

    hid_t dset  = H5Dcreate(group, name.c_str(), H5T_NATIVE_INT, filespace,
                      H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    
    H5Dwrite(dset, H5T_NATIVE_INT, memspace, filespace, dxpl_id, data);
    H5Dclose(dset);
    H5Sclose(memspace);
    H5Sclose(filespace);
    delete[] data;

}


void Writer::writeAttributeDbl(hid_t fileID, hid_t group, hid_t dxpl_id,
                            string name , double val){
    int rank ;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    
    int coresNum;
    MPI_Comm_size(MPI_COMM_WORLD, &coresNum);
    
    hsize_t offset, nLoc, nGlob;
    offset = rank;
    nLoc = 1;
    nGlob = coresNum;
    
    double * data = new double[1];
    data[0] = val;
    hid_t memspace  = H5Screate_simple(1, &nLoc , NULL);
    hid_t filespace = H5Screate_simple(1, &nGlob, NULL);
    H5Sselect_hyperslab(filespace, H5S_SELECT_SET, &offset, NULL, &nLoc, NULL);
    
    hid_t dset  = H5Dcreate(group, name.c_str(), H5T_NATIVE_DOUBLE, filespace,
                            H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    
    H5Dwrite(dset, H5T_NATIVE_DOUBLE, memspace, filespace, dxpl_id, data);
    H5Dclose(dset);
    H5Sclose(memspace);
    H5Sclose(filespace);
    delete[] data;
    
}

void Writer::writeAllForRestart(){
    
    int rank ;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    
    auto start_time = high_resolution_clock::now();
    
    logger->writeMsg("[Writer] writing restart ...", DEBUG);
    int subdomainXSize = loader->resolution[0];
    int subdomainYSize = loader->resolution[1];
    int subdomainZSize = loader->resolution[2];
    
    int idx;
    
    int totalNodeNumG1 = (subdomainXSize+1)*(subdomainYSize+1)*(subdomainZSize+1);
    int totalNodeNumG2 = (subdomainXSize+2)*(subdomainYSize+2)*(subdomainZSize+2);
    
    MPI_Info info = MPI_INFO_NULL;
    
    int coresNum;
    MPI_Comm_size(MPI_COMM_WORLD, &coresNum);
    
    int* gridSizesPerCore = new int[coresNum*3];
    int* particlesPerCore = new int[coresNum];
    
    MPI_Allgather(loader->resolution, 3, MPI_INT, gridSizesPerCore, 3, MPI_INT, MPI_COMM_WORLD);
    
    int totalPrtclNumber = pusher->getTotalParticleNumber();
    
    MPI_Allgather(&totalPrtclNumber, 1, MPI_INT, particlesPerCore, 1, MPI_INT, MPI_COMM_WORLD);
    
    int g2offset = 0;
    int g1offset = 0;
    int partoffset = 0;
    int g2glob = 0;
    int g1glob = 0;
    int partglob = 0;
    for ( int s=0; s<coresNum; s++){
        
        int lx = gridSizesPerCore[3*s+0];
        int ly = gridSizesPerCore[3*s+1];
        int lz = gridSizesPerCore[3*s+2];
        
        g1glob += (lx+1)*(ly+1)*(lz+1);
        g2glob += (lx+2)*(ly+2)*(lz+2);
        
        partglob += particlesPerCore[s];
        
        if(s<rank){
            g2offset += (lx+2)*(ly+2)*(lz+2);
            g1offset += (lx+1)*(ly+1)*(lz+1);
            
            partoffset += particlesPerCore[s];
        }
    }
    
    delete[] gridSizesPerCore;
    delete[] particlesPerCore;
    
    hid_t access = H5Pcreate(H5P_FILE_ACCESS);
    H5Pset_fapl_mpio(access, MPI_COMM_WORLD, info);
    
    string fileName = outputDir + fileNamePattern +  "restart.h5";
    hid_t fileID = H5Fcreate(fileName.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, access);
    
    const string groupname = "/vars";
    hid_t group   = H5Gcreate(fileID, groupname.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    
    hid_t dxpl_id = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(dxpl_id, H5FD_MPIO_COLLECTIVE);
    
    writeAttributeInt(fileID, group, dxpl_id, "g2_num", totalNodeNumG2);
    writeAttributeInt(fileID, group, dxpl_id, "g2_offset", g2offset);
    
    
    int totVarsOnG2 = gridMgr->getVarsNumOnG2();
    for ( int varN=0; varN<totVarsOnG2; varN++){
        VectorVar** vars = gridMgr->getVectorVariableOnG2(varN);
        int varSize = vars[0]->getSize();
        for ( int dir=0; dir<varSize; dir++){
            double* field = new double[totalNodeNumG2];
            for (idx = 0; idx < totalNodeNumG2; idx++) {
                field[idx] = vars[idx]->getValue()[dir];
            }
            string varName = "g2_"+to_string(varN)+"_"+to_string(dir);
            
            writeParallelWithOffset(fileID, group, dxpl_id, varName, field,
                                    totalNodeNumG2, g2glob, g2offset);
            delete[] field;
        }
    }

    writeAttributeInt(fileID, group, dxpl_id, "g1_num", totalNodeNumG1);
    writeAttributeInt(fileID, group, dxpl_id, "g1_offset", g1offset);
    
    int totVarsOnG1 = gridMgr->getVarsNumOnG1();
    for ( int varN=0; varN<totVarsOnG1; varN++){
        VectorVar** vars = gridMgr->getVectorVariableOnG1(varN);
        int varSize = vars[0]->getSize();
        for ( int dir=0; dir<varSize; dir++){
            double* field = new double[totalNodeNumG1];
            for (idx = 0; idx < totalNodeNumG1; idx++) {
                field[idx] = vars[idx]->getValue()[dir];
            }
            string varName = "g1_"+to_string(varN)+"_"+to_string(dir);
            writeParallelWithOffset(fileID, group, dxpl_id, varName, field,
                                    totalNodeNumG1, g1glob, g1offset);
            delete[] field;
        }
    }


    Particle** particles = pusher->getParticles();
    totalPrtclNumber = pusher->getTotalParticleNumber();
    
    writeAttributeInt(fileID, group, dxpl_id, "parts_num", totalPrtclNumber);
    writeAttributeInt(fileID, group, dxpl_id, "parts_offset", partoffset);
    
    int numOfSpecies = loader->getNumberOfSpecies();
    
    for(int spn = 0; spn<numOfSpecies;spn++){
        double w = pusher->getParticleWeight4Type(spn);
        writeAttributeDbl(fileID, group, dxpl_id, ("weight_"+to_string(spn)).c_str(), w);
    }
    
    double* particles2save = new double[PARTICLES_SIZE*totalPrtclNumber];
    int ind = 0;
    for ( idx=0; idx<totalPrtclNumber; idx++){
        particles[idx]->serialize(particles2save, PARTICLES_SIZE*idx);
    }
    
    writeParallelWithOffset(fileID, group, dxpl_id, "parts", particles2save,
                            PARTICLES_SIZE*totalPrtclNumber,
                            PARTICLES_SIZE*partglob,
                            PARTICLES_SIZE*partoffset);

    
    H5Fflush(fileID, H5F_SCOPE_GLOBAL);
    
    delete[] particles2save;
    
    H5Gclose(group);
    H5Pclose(dxpl_id);
   
    H5Fclose(fileID);
    
    
    auto end_time = high_resolution_clock::now();
    string msg ="[Writer] writing restart duration = "
    +to_string(duration_cast<milliseconds>(end_time - start_time).count())+" ms";
    logger->writeMsg(msg.c_str(), DEBUG);
}




