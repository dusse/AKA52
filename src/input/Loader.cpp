#include "Loader.hpp"

using namespace std;

Loader::Loader(){
    
}

const char  *INIT_CLASS_NAME = "Initializer";
const string  BRACKETS ="()";
const string  BRACKETS_3DOUBLE = "(ddd)";
const string  GET = "get";
const string  GET_TIMESTEP = "getTimestep";
const string  GET_MAX_TIMESTEPS_NUM = "getMaxTimestepsNum";
const string  GET_NUM_OF_SPECIES = "getNumOfSpecies";
const string  GET_TIMESTEP_WRITE = "getOutputTimestep";
const string  GET_OUTPUT_DIR = "getOutputDir";
const string  GET_FILENAME_TEMPLATE = "getOutputFilenameTemplate";
const string  GET_DENSITY = "getDensity";
const string  GET_NUM_OF_PARTICLES = "getParticlesPerCellNumber";
const string  GET_BFIELD   = "getBfield";
const string  GET_ELEPRES  = "getElectronPressure";
const string  GET_PRESSURE_SMOOTH_STRIDE = "getElectronPressureSmoothingStride";
const string  GET_VELOCITY = "getVelocity";
const string  GET_MASS   = "getMass";
const string  GET_CHARGE = "getCharge";
const string  SPIECIES   = "4species";
const string  GET_MPI_DOMAIN_NUM = "mpiDomainNum";

const string  GET_MIN_DENS_4_PPC = "getMinimumDens2ResolvePPC";

const string  GET_FIELD_BC_TYPE = "getFieldBCType";
const string  GET_PARTCL_BC_TYPE = "getParticleBCType";

const string  GET_RUN_TYPE = "getRunType";
const string  GET_INPUT_FILE = "getInputFile";

const string  GET_HYPERVISCOSITY = "getHyperviscosity";
const string  GET_ELECTRON_MASS  = "getElectronMass";
const string  GET_RELAX_FACTOR  = "getRelaxFactor";


const string  NUMBER_OF_LASER_SPOTS = "getNumberOfLaserSpots";
const string  GET_CENTERS  = "getLaserSpotCenter";
const string  GET_RADIUS   = "getLaserSpotRadius";

const string  PARTICLE_TYPE2HEAT = "getParticleType2Heat";
const string  PARTICLE_TEMP2LOAD = "getParticleTemp2Load";
const string  PARTICLE_DENS2KEEP = "getParticleDensity2sustain";
const string  PARTICLE_TEMP2KEEP = "getParticleTemperature2sustain";

const string  dirs[] = {"X", "Y", "Z"};

void Loader::load()
{
        MPI_Comm com;
        int rank ;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        int cs[3];
        int ss[3];

        Py_Initialize();
        this->pInstance = getPythonClassInstance(INIT_CLASS_NAME);


        for (int n=0; n<3; n++){
             string domainNum =GET+dirs[n]+GET_MPI_DOMAIN_NUM;
             this->mpiDomains[n] = (int) PyInt_AsLong(PyObject_CallMethod(pInstance,
                                                                           strdup(domainNum.c_str()),
                                                                           strdup(BRACKETS.c_str())));
            string bcType =GET_FIELD_BC_TYPE+dirs[n];
            this->BCtype[n] = (int) PyInt_AsLong(PyObject_CallMethod(pInstance,
                                                                         strdup(bcType.c_str()),
                                                                         strdup(BRACKETS.c_str())));// 1 - periodic for MPI
            string partclbcType =GET_PARTCL_BC_TYPE+dirs[n];
            this->partclBCtype[n] = (int) PyInt_AsLong(PyObject_CallMethod(pInstance,
                                                                 strdup(partclbcType.c_str()),
                                                                 strdup(BRACKETS.c_str())));// 1 - periodic, 0 - outflow
            switch (partclBCtype[n]) {
                case 0:
                    this->partclBCtype[n] = OUTFLOW_BC;
                    break;
                case 1:
                    this->partclBCtype[n] = PERIODIC_BC;
                    break;
                case 2:
                    this->partclBCtype[n] = REFLECT_BC;
                    break;
                    
                default:
                    throw runtime_error("no such particle BC!");
            }
            

            
        }
    
        string importantmsg0 = "[Loader] [BC] partclBCtype[0] = "+to_string( partclBCtype[0] )
                                    +" partclBCtype[1] = "+to_string(partclBCtype[1] )
                                    +" partclBCtype[2] = "+to_string(partclBCtype[2] );
    
        logger.writeMsg(importantmsg0.c_str(), INFO);
    
        string importantmsg00 = "[Loader] [BC] BCtype[0] = "+to_string( BCtype[0] )
                                +" BCtype[1] = "+to_string(BCtype[1] )
                                +" BCtype[2] = "+to_string(BCtype[2] );
    
        logger.writeMsg(importantmsg00.c_str(), INFO);

    
    
  
    
    
        this->runType = (int) PyInt_AsLong(PyObject_CallMethod(pInstance,
                                                           strdup(GET_RUN_TYPE.c_str()),
                                                           strdup(BRACKETS.c_str())));;
    
        this->inputfile = PyString_AS_STRING(PyObject_CallMethod(pInstance,
                                                             strdup(GET_INPUT_FILE.c_str()),
                                                             strdup(BRACKETS.c_str())));;
    
        MPI_Cart_create(MPI_COMM_WORLD, 3, mpiDomains, BCtype, 1, &com);
        MPI_Cart_coords(com, rank, 3, cs);
    
    string importantmsg = "[Loader] [MPI] mpiDomain[0] = "+to_string(mpiDomains[0])
    +" mpiDomain[1] = "+to_string(mpiDomains[1])
    +" mpiDomain[2] = "+to_string(mpiDomains[2]);
    
    logger.writeMsg(importantmsg.c_str(), INFO);
    
        for(int i = 0; i < 3; i++){
            this->mpiCoords[i] = cs[i];
            this->BCtype[i] = BCtype[i] == 1 ? PERIODIC : IDEAL;
        }
        int neighborRank;

        for(int i = 0; i < 27; i++){
            this->neighbors2Send.push_back(MPI_PROC_NULL);
            this->neighbors2Recv.push_back(MPI_PROC_NULL);
        }
    
        for (int a = -1; a <= +1; a++){
            for (int b = -1; b <= +1; b++){
                for (int c = -1; c <= +1; c++){
                    
                    //-1 -> left, bottom, back
                    // 0 -> same position
                    //+1 -> right, top, front
                    
                    /* __ guess the coordinate of the neighbor __ */
                    ss[0] = cs[0]+a;
                    ss[1] = cs[1]+b;
                    ss[2] = cs[2]+c;
                    
                    
                    /* __ if coordinate out of range : no neighbor __ */
                    if ( (BCtype[0] != PERIODIC
                              && (ss[0] < 0 || ss[0] >= mpiDomains[0]))
                          || ( BCtype[1] != PERIODIC
                              && (ss[1] < 0 || ss[1] >= mpiDomains[1]))
                          || ( BCtype[2] != PERIODIC
                              && (ss[2] < 0 || ss[2] >= mpiDomains[2]))){
                          
                          neighborRank = MPI_PROC_NULL;
                              
                        //  For a process group with cartesian structure,
                        //  the function MPI_CART_RANK translates the logical
                        //  process coordinates to process ranks as they are used
                        //  by the point-to-point routines.
                        //  For dimension i with periods(i) = true, if the coordinate,
                        //  coords(i), is out of range, that is, coords(i) < 0 or coords(i) >= dims(i),
                        //  it is shifted back to the interval 0 <= coords(i) < dims(i) automatically.
                        //
                        //  Out-of-range coordinates are erroneous for non-periodic dimensions.
                        //  Versions of MPICH before 1.2.2 returned MPI_PROC_NULL for the rank in this case.
                              
                    }else{
                          
                          MPI_Cart_rank(com, ss, &neighborRank);
                          
                    }

                    this->neighbors2Send[(1+c)+3*((1+b)+3*(1+a))] = neighborRank;
                    this->neighbors2Recv[(1-c)+3*((1-b)+3*(1-a))] = neighborRank;
                    
                    }
                }
            }


        double boxSizePerDomain[3];
        for (int n=0; n<3; n++){
        
            string res = GET + dirs[n] + "resolution";
            string right_str= GET + dirs[n] + "right";
            this->boxCoordinates[n][0] = 0.0; //origin point (0,0,0) is left lower corner i=0 j=0 k=0
            this->boxCoordinates[n][1] = PyFloat_AsDouble(PyObject_CallMethod(pInstance,
                                                                              strdup(right_str.c_str()),
                                                                              strdup(BRACKETS.c_str())));
            this->boxSizes[n] = boxCoordinates[n][1];
            //  length of the box in normalized units
            
            if(boxSizes[n] <= 0.0){
                throw runtime_error("box has zero length!");
            }
            
            //  number of pixel per box
            this->totPixelsPerBoxSide[n] = (int) PyInt_AsLong(PyObject_CallMethod(pInstance,
                                                                                  strdup(res.c_str()),
                                                                                  strdup(BRACKETS.c_str())));
            
            if (totPixelsPerBoxSide[n] == 1){
                // length per pixel equals to total size
                this->spatialSteps[n]=boxSizes[n];
            }else{
                this->spatialSteps[n]=boxSizes[n]/(double)(totPixelsPerBoxSide[n]);
            }
            
            double approxRes = int(totPixelsPerBoxSide[n]/double(mpiDomains[n]));
            
            this->resolution[n] = int(approxRes);
            
            int delta = totPixelsPerBoxSide[n] - resolution[n]*mpiDomains[n];
            
            this->offsetInPixels[n] = 0;
            
            if(delta == 0.0){
                //  length of the domain in normalized units
                boxSizePerDomain[n]  = resolution[n]*spatialSteps[n];
                this->offsetInPixels[n]    = resolution[n]*cs[n];// global offset in pixels
                
            }else if(cs[n] < delta){
                
                this->resolution[n] += 1;
                
                boxSizePerDomain[n]  = resolution[n]*spatialSteps[n];
                this->offsetInPixels[n]    = resolution[n]*cs[n];
                
            } else {
                boxSizePerDomain[n] = resolution[n]*spatialSteps[n];
                this->offsetInPixels[n] = (resolution[n]+1)*delta + resolution[n]*(cs[n]-delta);
            }
            
            this->boxCoordinates[n][0] = offsetInPixels[n]*spatialSteps[n];
            this->boxCoordinates[n][1] = boxCoordinates[n][0]+boxSizePerDomain[n];
            
            string msgs = "[Loader] [MPI] rank = "+to_string(rank)
            +"\n"+string( 10, ' ' )+" resolution = "+to_string(resolution[n])
            +"\n"+string( 10, ' ' )+" offsetInPixels = "+to_string(offsetInPixels[n])
            +"\n"+string( 10, ' ' )+" boxCoordinates[n][0] = "+to_string(boxCoordinates[n][0])
            +"\n"+string( 10, ' ' )+" boxCoordinates[n][1] = "+to_string(boxCoordinates[n][1]);
            logger.writeMsg(msgs.c_str(), INFO);
        
    }
    if(boxSizes[0] != 1 && boxSizes[1] != 1 && boxSizes[2] != 1 ){
        this->dim = 3;
    }else if(boxSizes[0] != 1 && boxSizes[1] != 1 ){
        this->dim = 2;
    }else if(boxSizes[0] != 1 ){
        this->dim = 1;
    }
    
    if(boxSizes[0] != 1 && boxSizes[1] != 1 && boxSizes[2] != 1 ){
        this->dim = 3;
    }
    
    this->timeStep           = PyFloat_AsDouble(PyObject_CallMethod(pInstance,
                                                          strdup(GET_TIMESTEP.c_str()),
                                                          strdup(BRACKETS.c_str())));
    
    this->numOfSpecies       = PyFloat_AsDouble(PyObject_CallMethod(pInstance,
                                                              strdup(GET_NUM_OF_SPECIES.c_str()),
                                                              strdup(BRACKETS.c_str())));
    
    this->ppc                = (int) PyFloat_AsDouble(PyObject_CallMethod(pInstance,
                                                           strdup(GET_NUM_OF_PARTICLES.c_str()),
                                                           strdup(BRACKETS.c_str())));
    
    this->minimumDens2ResolvePPC =  PyFloat_AsDouble(PyObject_CallMethod(pInstance,
                                                                              strdup(GET_MIN_DENS_4_PPC.c_str()),
                                                                              strdup(BRACKETS.c_str())));

    
    this->maxTimestepsNum    = PyFloat_AsDouble(PyObject_CallMethod(pInstance,
                                                                 strdup(GET_MAX_TIMESTEPS_NUM.c_str()),
                                                                 strdup(BRACKETS.c_str())));
    
    this->timestepsNum2Write = PyFloat_AsDouble(PyObject_CallMethod(pInstance,
                                                                    strdup(GET_TIMESTEP_WRITE.c_str()),
                                                                    strdup(BRACKETS.c_str())));
    
    this->outputDir          = PyString_AS_STRING(PyObject_CallMethod(pInstance,
                                                           strdup(GET_OUTPUT_DIR.c_str()),
                                                           strdup(BRACKETS.c_str())));
    
    this->fileNameTemplate   = PyString_AS_STRING(PyObject_CallMethod(pInstance,
                                                                  strdup(GET_FILENAME_TEMPLATE.c_str()),
                                                                  strdup(BRACKETS.c_str())));
    
    this->hyperviscosity     = PyFloat_AsDouble(PyObject_CallMethod(pInstance,
                                                                  strdup(GET_HYPERVISCOSITY.c_str()),
                                                                  strdup(BRACKETS.c_str())));
    
    this->electronmass       = PyFloat_AsDouble(PyObject_CallMethod(pInstance,
                                                                    strdup(GET_ELECTRON_MASS.c_str()),
                                                                    strdup(BRACKETS.c_str())));
    
    this->smoothStride       = (int) PyFloat_AsDouble(PyObject_CallMethod(pInstance,
                                                                    strdup(GET_PRESSURE_SMOOTH_STRIDE.c_str()),
                                                                    strdup(BRACKETS.c_str())));
    
    this->relaxFactor        = PyFloat_AsDouble(PyObject_CallMethod(pInstance,
                                                                    strdup(GET_RELAX_FACTOR.c_str()),
                                                                    strdup(BRACKETS.c_str())));
    
    
    this->numOfSpots         = PyFloat_AsDouble(PyObject_CallMethod(pInstance,
                                                            strdup(NUMBER_OF_LASER_SPOTS.c_str()),
                                                            strdup(BRACKETS.c_str())));
    
    if(numOfSpots>0){
        
        this->prtclType2Heat = PyFloat_AsDouble(PyObject_CallMethod(pInstance,
                                                                   strdup(PARTICLE_TYPE2HEAT.c_str()),
                                                                   strdup(BRACKETS.c_str())))-1;
        // -1 because in input count starts in human manner from 1
        
        this->prtclTemp2Load = PyFloat_AsDouble(PyObject_CallMethod(pInstance,
                                                                    strdup(PARTICLE_TEMP2LOAD.c_str()),
                                                                    strdup(BRACKETS.c_str())));
        
        
        this->prtclDens2Keep = PyFloat_AsDouble(PyObject_CallMethod(pInstance,
                                                                   strdup(PARTICLE_DENS2KEEP.c_str()),
                                                                   strdup(BRACKETS.c_str())));
        
        this->prtclTemp2Keep = PyFloat_AsDouble(PyObject_CallMethod(pInstance,
                                                                   strdup(PARTICLE_TEMP2KEEP.c_str()),
                                                                   strdup(BRACKETS.c_str())));
        
        
        centersOfSpots = new double[numOfSpots*3*sizeof(double)];
        spotsRadius    = new double[numOfSpots*sizeof(double)];
        
        for (int n=0; n<numOfSpots; n++){
            
            for(int dim=0; dim<3; dim++){
                string varName =GET_CENTERS+dirs[dim]+"position"+to_string(n+1);
                centersOfSpots[3*n+dim] = PyFloat_AsDouble(PyObject_CallMethod(pInstance,
                                                                               strdup(varName.c_str()),
                                                                               strdup(BRACKETS.c_str())));
            }
            
            string varName =GET_RADIUS+to_string(n+1);
            spotsRadius[n] = PyFloat_AsDouble(PyObject_CallMethod(pInstance,
                                                                  strdup(varName.c_str()),
                                                                  strdup(BRACKETS.c_str())));
        }
    }

    
    if(rank == 0){

        string msg = runType == 0 ? "scratch" : "restart";
        msg = "[Loader] [COMMON] run from "+msg+" "+to_string(dim)+"D simulation";
        logger.writeMsg(msg.c_str(), INFO);
        
        msg = "[Loader] [COMMON]  Box size:";
        logger.writeMsg(msg.c_str(), INFO);

        for (int n=0; n<3; n++){
                   char buf[100];
                   snprintf(buf, sizeof(buf),
                            "[Loader] [COMMON] [%s] [%1.5f, %1.5f] res = %3i => step  = %1.6f",
                   dirs[n].c_str(), boxCoordinates[n][0],boxSizes[n],
                   totPixelsPerBoxSide[n],spatialSteps[n]);
            
            
            logger.writeMsg(buf, INFO);
        }
        msg = "[Loader] [COMMON] timeStep = "+to_string(timeStep);
        logger.writeMsg(msg.c_str(), INFO);
        
        
        msg = "[Loader] [COMMON] Max timeStep number = "+to_string(maxTimestepsNum);
        logger.writeMsg(msg.c_str(), INFO);
        msg = "[Loader] [COMMON] timeStep number to write to file = "+to_string(timestepsNum2Write);
        logger.writeMsg(msg.c_str(), INFO);
        msg = "[Loader] [COMMON] outputDir = "+outputDir;
        logger.writeMsg(msg.c_str(), INFO);
        msg = "[Loader] [COMMON] filename template = "+fileNameTemplate;
        logger.writeMsg(msg.c_str(), INFO);
        msg = "[Loader] [COMMON] number of species = "+to_string(numOfSpecies);
        logger.writeMsg(msg.c_str(), INFO);
        msg = "[Loader] [COMMON] ppc = "+to_string(ppc);
        logger.writeMsg(msg.c_str(), INFO);
        msg = "[Loader] [COMMON] minimum density resolved by ppc number = "+to_string(minimumDens2ResolvePPC);
        logger.writeMsg(msg.c_str(), INFO);
        
        if(numOfSpots>0){
            msg = "[Loader] [LASER] numOfSpots = "+to_string(numOfSpots)
                    +"; prtclType2Heat = "+to_string(prtclType2Heat+1)
                    +"; prtclDens2Keep = "+to_string(prtclDens2Keep)
                    +"; prtclTemp2Keep = "+to_string(prtclTemp2Keep);
            logger.writeMsg(msg.c_str(), INFO);
        
            for (int n=0; n<numOfSpots; n++){
                char buf[100];
                snprintf(buf, sizeof(buf),"[Loader] [LASER] spot [%s] [%s, %s, %s] R = %f",
                       to_string(n).c_str(),
                       to_string(centersOfSpots[3*n+0]).c_str(),
                       to_string(centersOfSpots[3*n+1]).c_str(),
                       to_string(centersOfSpots[3*n+2]).c_str(),
                       spotsRadius[n]);
                logger.writeMsg(buf, INFO);
            }
        }
        msg = "[Loader] [OHM's LAW]: hyperviscosity = "+to_string(hyperviscosity);
        logger.writeMsg(msg.c_str(), INFO);
        
        msg = "[Loader] [PRESSURE]: electron mass = "+to_string(electronmass);
        logger.writeMsg(msg.c_str(), INFO);

        msg = "[Loader] [PRESSURE]: relaxation factor = "+to_string(relaxFactor);
        logger.writeMsg(msg.c_str(), INFO);
        
        msg = "[Loader] [PRESSURE]: smoothing stride = "+to_string(smoothStride);
        logger.writeMsg(msg.c_str(), INFO);
        
        logger.writeMsg("[Loader] initialization...OK!", DEBUG);
    }

}

int Loader::getNumberOfSpecies(){
    return numOfSpecies;
}

double Loader::getTimeStep(){
    return timeStep;
}

int Loader::getMaxTimestepsNum(){
    return maxTimestepsNum;
}

int Loader::getTimestepsNum2Write(){
    return timestepsNum2Write;
}

string Loader::getOutputDir() const {
    return outputDir;
}

string Loader::getFilenameTemplate(){
    return fileNameTemplate;
}

vector<double> Loader::getVelocity(double x, double y, double z,
                                   int spieceiesType){
    PyObject *pValue;
    vector<double> velocity;
    
    for (int n=0; n<3; n++){
        string varName =GET_VELOCITY+dirs[n]+SPIECIES+to_string(spieceiesType+1);

        pValue = PyObject_CallMethod(pInstance, strdup(varName.c_str()),
                                     strdup(BRACKETS_3DOUBLE.c_str()),x,y,z);
        velocity.push_back(PyFloat_AsDouble(pValue));
    }
    return velocity;
}




vector<double> Loader::getBfield(double x,double y,double z){
    PyObject *pValue;
    vector<double> bField;
    
    for (int n=0; n<3; n++){
        string varName =GET_BFIELD+dirs[n];
        
        pValue = PyObject_CallMethod(pInstance, strdup(varName.c_str()),
                                     strdup(BRACKETS_3DOUBLE.c_str()),x,y,z);
        bField.push_back(PyFloat_AsDouble(pValue));
    }
    return bField;
}

double Loader::getDensity(double x,double y,double z, int spieceiesType){
    string varName = GET_DENSITY+SPIECIES+to_string(spieceiesType+1);
    return PyFloat_AsDouble(PyObject_CallMethod(pInstance,
                                                strdup(varName.c_str()),
                                                strdup(BRACKETS_3DOUBLE.c_str()),
                                                x,y,z));
}


double Loader::getMass4spieceies(int spieceiesType){
    string varName = GET_MASS+SPIECIES+to_string(spieceiesType+1);
    return PyFloat_AsDouble(PyObject_CallMethod(pInstance,
                                                strdup(varName.c_str()),
                                                strdup(BRACKETS.c_str())));
}

double Loader::getCharge4spieceies(int spieceiesType){
    string varName = GET_CHARGE+SPIECIES+to_string(spieceiesType+1);
    return PyFloat_AsDouble(PyObject_CallMethod(pInstance,
                                                strdup(varName.c_str()),
                                                strdup(BRACKETS.c_str())));
}


double Loader::getElectronPressure(double x,double y,double z){
    string varName = GET_ELEPRES;
    return PyFloat_AsDouble(PyObject_CallMethod(pInstance,
                                                strdup(varName.c_str()),
                                                strdup(BRACKETS_3DOUBLE.c_str()),
                                                x,y,z));
}


double Loader::getParticlesPerCellNumber(){
    return ppc;
}

PyObject * Loader::getPythonClassInstance(string className){
    PyObject  *pName, *pModule, *pDict, *pClass, *pInstance;
    string msg = "[Loader] Start to instantiate the Python class " + className;
    logger.writeMsg(msg.c_str(), DEBUG);
    pName = PyString_FromString(className.c_str());
    pModule = PyImport_Import(pName);
    pDict = PyModule_GetDict(pModule);
    pClass = PyDict_GetItemString(pDict, className.c_str());
    if (PyCallable_Check(pClass))
    {
        pInstance = PyObject_CallObject(pClass, NULL);
    }
    else
    {
        logger.writeMsg("[Loader] Cannot instantiate the Python class", CRITICAL);
        pInstance = nullptr;
    }
    logger.writeMsg("[Loader] finish to instantiate the Python class ", DEBUG);
    
    return pInstance;
}

Loader::~Loader()
{
    Py_Finalize();
    logger.writeMsg("[Loader] FINALIZE...OK!", DEBUG);
}
