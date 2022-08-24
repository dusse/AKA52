#include "Loader.hpp"

using namespace std;

Loader::Loader(){
    
}

/*#################################################################################################*/

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
const string  GET_BFIELD   = "getBfield";
const string  GET_ELEPRESXX  = "getElectronPressureXX";
const string  GET_ELEPRESYY  = "getElectronPressureYY";
const string  GET_ELEPRESZZ  = "getElectronPressureZZ";
const string  GET_PRESSURE_SMOOTH_STRIDE = "getElectronPressureSmoothingStride";
const string  GET_VELOCITY = "getVelocity";
const string  GET_FLUID_VELOCITY = "getFluidVelocity";
const string  INJECTED_PARTICLES = "4InjectedParticles";
const string  GET_DFTYPE = "getDFtype";
const string  GET_DFTYPE4INJECTEDPARTICLES = "getDFtype4InjectedParticles";
const string  GET_MASS   = "getMass";
const string  GET_CHARGE = "getCharge";
const string  GET_IFPARTICLETYPEISFROZEN = "getIfParticleTypeIsFrozen";
const string  GET_PPC = "getPPC";
const string  SPECIES   = "4species";
const string  GET_PPC4LOADED_PARTICLES = "getPPC4loadedParticles";
const string  GET_MPI_DOMAIN_NUM = "mpiDomainNum";

const string  GET_MIN_DENS_4_PPC = "getMinimumDens2ResolvePPC";

const string  GET_FIELD_BC_TYPE = "getFieldBCType";
const string  GET_PARTCL_BC_TYPE = "getParticleBCType";

const string  GET_DAMPING_BOUNDARY_WIDTH  = "getDampingBoundaryWidth";

const string  LEFT  = "left";
const string  RIGHT  = "right";

const string  GET_RUN_TYPE = "getRunType";
const string  GET_INPUT_FILE = "getInputFile";

const string  GET_RESISTIVITY = "getResistivity";
const string  GET_BFIELD_LIMIT    = "getBfieldLimit";

const string  GET_DEFAULT_COLOUMB_LOG = "getDefaultColoumbLogarithm";
const string  GET_IONION_COLLISION_FREQUENCY_FACTOR = "getIonIonCollisionFrequencyFactor";

const string  GET_ELECTRON_MASS  = "getElectronMass";
const string  GET_RELAX_FACTOR   = "getRelaxFactor";

const string  IF2USE_ISOTHERMAL_CLOSURE = "getIfWeUseIsothermalClosure";
const string  GET_FIXED_ELECTRON_TEMPERATURE = "getElectronTemperature4IsothermalClosure";

const string  GET_CELL_BREAKDOWN_EFIELD_FACTOR  = "getCellBreakdownEfieldFactor";
const string  GET_CRITICAL_PRESSURE   = "getCriticalPressure";

const string  NUMBER_OF_LASER_SPOTS = "getNumberOfLaserSpots";

const string  PARTICLE_TYPE2LOAD = "getParticleType2Load";

const string  GET_LASER_PULSE_DURATION = "getLaserPulseDuration";


const string  GET_TARGET_ION_DENSITY2SUSTAIN = "getTargetIonDensity2sustain";
const string  GET_ELECTRON_PRESSURE2SUSTAIN  = "getElectronPressure2sustain";

const string  dirs[] = {"X", "Y", "Z"};


/*#################################################################################################*/

PyObject * Loader::getPythonClassInstance(string className){
    PyObject  *pName, *pClass, *pModule, *pDict;
    string msg = "[Loader] Start to instantiate the Python class " + className;
    logger.writeMsg(msg.c_str(), DEBUG);
    
    string inputFilePath = string(getenv("INPUTFILEPATH"))+ INIT_CLASS_NAME + ".py";
    string introMsg = "[Loader] used input file is " + inputFilePath;
    logger.writeMsg(introMsg.c_str(), INFO);

#if PY_MAJOR_VERSION >= 3
    pName = PyUnicode_FromString(className.c_str());
#else
    pName = PyString_FromString(className.c_str());
#endif
    
    pModule = PyImport_Import(pName);
    
    if( pModule == NULL ){
        logger.writeMsg("*****************************************************", CRITICAL);
        logger.writeMsg("****                                             ****", CRITICAL);
        logger.writeMsg("****  STOP SIMULATION!!!    INPUT FILE PROBLEM   ****", CRITICAL);
        logger.writeMsg("****                                             ****", CRITICAL);
        logger.writeMsg("****   try debug mode 'make -DLOG'               ****", CRITICAL);
        logger.writeMsg("****   check indents in python input file        ****", CRITICAL);
        logger.writeMsg("****   check python enviroment and imports       ****", CRITICAL);
        logger.writeMsg("*****************************************************", CRITICAL);
        exit(-1);
    }
    
    pDict = PyModule_GetDict(pModule);
    pClass = PyDict_GetItemString(pDict, className.c_str());
    
    if( PyCallable_Check(pClass) ){
        pInstance = PyObject_CallObject(pClass, NULL);
    }else{
        logger.writeMsg("[Loader] Cannot instantiate the Python class", CRITICAL);
        pInstance = nullptr;
    }
    
    logger.writeMsg("[Loader] finish to instantiate the Python class ", DEBUG);
    
    return pInstance;
}



double Loader::callPyFloatFunction( PyObject* instance,
                                    const string funcName,
                                    const string brackets){

    if (checkMethodExistence(funcName, brackets) == METHOD_OK){
        return  PyFloat_AsDouble(getPyMethod(instance,funcName,brackets));
    }else{
        return 0.0;
    }
}

double Loader::callPyFloatFunctionWith3args( PyObject* instance,
                                            const string funcName,
                                            const string brackets,
                                            double x,double y,double z){

    if (checkMethodExistence(funcName, brackets) == METHOD_OK){
        return  PyFloat_AsDouble(PyObject_CallMethod(instance,strdup(funcName.c_str()),
                                                strdup(brackets.c_str()),x,y,z));
    }else{
        return 0.0;
    } 
}

long Loader::callPyLongFunction( PyObject* instance,
                                const string funcName,
                                const string brackets){

    if (checkMethodExistence(funcName, brackets) == METHOD_OK){
	#if PY_MAJOR_VERSION >= 3
 	return PyLong_AS_LONG(getPyMethod(instance,funcName,brackets));
	#else
 	return PyInt_AsLong(getPyMethod(instance,funcName,brackets));
	#endif
        
    }else{
 	return 0.0;
    }    
}


string Loader::callPyStringFunction( PyObject* instance,
                                    const string funcName,
                                    const string brackets){
    
#if PY_MAJOR_VERSION >= 3
    PyObject* str = PyUnicode_AsUTF8String(getPyMethod(instance,funcName,brackets));
    return PyBytes_AsString(str);
#else
    return PyString_AS_STRING(getPyMethod(instance,funcName,brackets));
#endif
    
}


PyObject* Loader::getPyMethod(PyObject* instance,
                               const string funcName,
                               const string brackets){
    return PyObject_CallMethod(instance, strdup(funcName.c_str()),
                               strdup(brackets.c_str()));
}


int Loader::checkMethodExistence(const string funcName, const string brackets){

    if( checkedMethods.count(funcName) ){
        return METHOD_OK;
    }

    if( failMethods.count(funcName) ){
        return METHOD_FAIL;
    }

    string inputFilePath = string(getenv("INPUTFILEPATH"))+ INIT_CLASS_NAME + ".py";

    ifstream inFile;
    inFile.open(inputFilePath);

    string line;

    while(getline(inFile, line)) {
        if( line.find(funcName, 0) != string::npos ){
	   if( brackets == BRACKETS_3DOUBLE && line.find(",", 0) == string::npos ){
		string msg = "[Loader] [METHOD SIGNATURE PROBLEM] method "+ funcName +" has inappropriate format, expected 3 doubles ";
                logger.writeMsg(msg.c_str(), CRITICAL); 
		throw runtime_error("METHOD SIGNATURE PROBLEM. CHECK INPUT FILE!!!");
	   }
	   checkedMethods.insert(funcName);
	   return METHOD_OK;
        }
    }

    string msg = "[Loader] [METHOD NOT FOUND] method "+ funcName +" is absent in Initializer.py. Used default value 0.";
    logger.writeMsg(msg.c_str(), INFO); 
    failMethods.insert(funcName);

    return METHOD_FAIL;
}

/*#################################################################################################*/

void Loader::initMPIcoordinatesOfDomains( int rank, int cs[3] ){
    
    int ss[3];
    MPI_Comm com;
    
    MPI_Cart_create(MPI_COMM_WORLD, 3, this->mpiDomains, this->BCtype, 1, &com);
    MPI_Cart_coords(com, rank, 3, cs);
    
    
    for( int i = 0; i < 3; i++ ){
        this->mpiCoords[i] = cs[i];
        this->BCtype[i] = BCtype[i] == 1 ? PERIODIC : DAMPING;
    }
    int neighborRank;
    
    for( int i = 0; i < 27; i++ ){
        this->neighbors2Send.push_back(MPI_PROC_NULL);
        this->neighbors2Recv.push_back(MPI_PROC_NULL);
    }
    
    //-1 -> left, bottom, back
    // 0 -> same position
    //+1 -> right, top, front
    
    for( int a = -1; a <= +1; a++ ){
        for( int b = -1; b <= +1; b++ ){
            for( int c = -1; c <= +1; c++ ){
                
                /* __ coordinate of the neighbor __ */
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
}


void Loader::load(){
    
        int rank ;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);

        Py_Initialize();
        this->pInstance = getPythonClassInstance(INIT_CLASS_NAME);


        for( int n = 0; n < 3; n++ ){
            
            string domainNum =GET+dirs[n]+GET_MPI_DOMAIN_NUM;
            this->mpiDomains[n] = (int) callPyLongFunction( pInstance, domainNum, BRACKETS );
            
            string bcType =GET_FIELD_BC_TYPE+dirs[n];// 1 - periodic for MPI, 0 - damping layer for B and Pe
            this->BCtype[n] = (int) callPyLongFunction( pInstance, bcType, BRACKETS);
            
            string partclbcType =GET_PARTCL_BC_TYPE+dirs[n];// 1 - periodic, 0 - outflow
            this->partclBCtype[n] = (int) callPyLongFunction( pInstance, partclbcType, BRACKETS );
            
            switch (partclBCtype[n]) {
                case 0:
                    this->partclBCtype[n] = OUTFLOW_BC;
                    break;
                case 1:
                    this->partclBCtype[n] = PERIODIC_BC;
                    break;
                default:
                    throw runtime_error("no such particle BC!");
            }
            
            string right_damping = GET_DAMPING_BOUNDARY_WIDTH + dirs[n] + RIGHT;
            string left_damping  = GET_DAMPING_BOUNDARY_WIDTH + dirs[n] + LEFT;
            
            PyObject* calldampingBoundaryWidthMethod = getPyMethod( pInstance, left_damping, BRACKETS );
            
            if( calldampingBoundaryWidthMethod != NULL ){
                this->dampingBoundaryWidth[n][0] = PyFloat_AsDouble(calldampingBoundaryWidthMethod);
            }
            
            calldampingBoundaryWidthMethod = getPyMethod( pInstance, right_damping, BRACKETS );
            
            if( calldampingBoundaryWidthMethod != NULL ){
                this->dampingBoundaryWidth[n][1] = PyFloat_AsDouble(calldampingBoundaryWidthMethod);
            }
            
        }
    
        string pbcmsg = "[Loader] [BC] partclBCtype[0] = "+to_string(partclBCtype[0])
                                    +" partclBCtype[1] = "+to_string(partclBCtype[1])
                                    +" partclBCtype[2] = "+to_string(partclBCtype[2]);
    
        logger.writeMsg(pbcmsg.c_str(), INFO);
    
        string fbcmsg = "[Loader] [BC] BCtype[0] = "+to_string(BCtype[0])
                                    +" BCtype[1] = "+to_string(BCtype[1])
                                    +" BCtype[2] = "+to_string(BCtype[2]);
    
        logger.writeMsg(fbcmsg.c_str(), INFO);

    
        if( this->BCtype[0] ==  DAMPING ){
            string dampingBCmsg = "[Loader] [BC - X] damping layer width: left = "
                                +to_string( this->dampingBoundaryWidth[0][0] )
                                +" / right = "+to_string(this->dampingBoundaryWidth[0][1] );
            logger.writeMsg(dampingBCmsg.c_str(), INFO);
        }else{
            this->dampingBoundaryWidth[0][0] = 0.0;
            this->dampingBoundaryWidth[0][1] = 0.0;
        }
    
        if( this->BCtype[1] ==  DAMPING ){
            string dampingBCmsg = "[Loader] [BC - Y] damping layer width: left = "
                                +to_string( this->dampingBoundaryWidth[1][0] )
                                +" / right = "+to_string(this->dampingBoundaryWidth[1][1] );
            logger.writeMsg(dampingBCmsg.c_str(), INFO);
        }else{
            this->dampingBoundaryWidth[1][0] = 0.0;
            this->dampingBoundaryWidth[1][1] = 0.0;
        }
    
        if( this->BCtype[2] ==  DAMPING ){
            string dampingBCmsg = "[Loader] [BC - Z] damping layer width: left = "
                                +to_string( this->dampingBoundaryWidth[2][0] )
                                +" / right = "+to_string(this->dampingBoundaryWidth[2][1] );
            logger.writeMsg(dampingBCmsg.c_str(), INFO);
        }else{
            this->dampingBoundaryWidth[2][0] = 0.0;
            this->dampingBoundaryWidth[2][1] = 0.0;
        }
  
    
    
        this->runType = (int) callPyLongFunction( pInstance, GET_RUN_TYPE, BRACKETS );
    
        this->inputfile =  callPyStringFunction( pInstance, GET_INPUT_FILE, BRACKETS );
    
        int nprocs;
        MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    
        string mpimsg = "[Loader] [MPI] mpiDomain[0] (X) = "+to_string(mpiDomains[0])
                                           +" mpiDomain[1] (Y) = "+to_string(mpiDomains[1])
                                           +" mpiDomain[2] (Z) = "+to_string(mpiDomains[2]);
    
        logger.writeMsg(mpimsg.c_str(), INFO);
    
        if( mpiDomains[0]*mpiDomains[1]*mpiDomains[2] != nprocs ){
            string errmsg = "[Loader] [MPI] number of cores run with mpirun command  = "+to_string(nprocs);
            logger.writeMsg(errmsg.c_str(), CRITICAL);
            throw runtime_error("!!!number of requested cores must be equal to number of distributed cores!!!");
        }
    
        int cs[3];
        initMPIcoordinatesOfDomains( rank, cs );

        double boxSizePerDomain[3];
    
        for( int n = 0; n < 3; n++ ){
        
            string res = GET + dirs[n] + "resolution";
            string right_str= GET + dirs[n] + RIGHT;
            this->boxCoordinates[n][0] = 0.0; //origin point (0,0,0) is left lower corner i=0 j=0 k=0
            this->boxCoordinates[n][1] = callPyFloatFunction( pInstance, right_str, BRACKETS );
            
            this->boxSizes[n] = boxCoordinates[n][1];
            //  length of the box in normalized units
            
            if( boxSizes[n] <= 0.0 ){
                throw runtime_error("box has zero length!");
            }
            
            //  number of pixel per box
            this->totPixelsPerBoxSide[n] = (int) callPyLongFunction( pInstance, res, BRACKETS);
            
            if( totPixelsPerBoxSide[n] == 1 ){
                // length per pixel equals to total size
                this->spatialSteps[n]=boxSizes[n];
            }else{
                this->spatialSteps[n]=boxSizes[n]/(double)(totPixelsPerBoxSide[n]);
            }
            
            double approxRes = int(totPixelsPerBoxSide[n]/double(mpiDomains[n]));
            
            this->resolution[n] = int(approxRes);
            
            int delta = totPixelsPerBoxSide[n] - resolution[n]*mpiDomains[n];
            
            this->offsetInPixels[n] = 0;
            
            if ( delta == 0.0 ){
                //  length of the domain in normalized units
                boxSizePerDomain[n]     = resolution[n]*spatialSteps[n];
                this->offsetInPixels[n] = resolution[n]*cs[n];// global offset in pixels
                
            } else if ( cs[n] < delta ){
                this->resolution[n] += 1;
                
                boxSizePerDomain[n]     = resolution[n]*spatialSteps[n];
                this->offsetInPixels[n] = resolution[n]*cs[n];
                
            } else {
                boxSizePerDomain[n]     = resolution[n]*spatialSteps[n];
                this->offsetInPixels[n] = (resolution[n]+1)*delta + resolution[n]*(cs[n]-delta);
            }
            
            this->boxCoordinates[n][0] = offsetInPixels[n]*spatialSteps[n];
            this->boxCoordinates[n][1] = boxCoordinates[n][0]+boxSizePerDomain[n];
            
            string msgs = "[Loader] [MPI] rank = "+to_string(rank)
            +"\n"+string( 10, ' ' )+" resolution = "+to_string(resolution[n])
            +"\n"+string( 10, ' ' )+" offsetInPixels = "+to_string(offsetInPixels[n])
            +"\n"+string( 10, ' ' )+" boxCoordinates["+to_string(n)+"][0] = "+to_string(boxCoordinates[n][0])
            +"\n"+string( 10, ' ' )+" boxCoordinates["+to_string(n)+"][1] = "+to_string(boxCoordinates[n][1]);
            logger.writeMsg(msgs.c_str(), INFO);

    }
    
    // needed for BoundaryManager.cpp
    if( boxSizes[0] != 1 && boxSizes[1] != 1 && boxSizes[2] != 1 ){
        this->dim = 3;
    }else if( boxSizes[0] != 1 && boxSizes[1] != 1 ){
        this->dim = 2;
    }else{
        this->dim = 1;
    }
    
    
    this->timeStep              =  callPyFloatFunction( pInstance, GET_TIMESTEP, BRACKETS );
    
    this->numOfSpecies          =  callPyFloatFunction( pInstance, GET_NUM_OF_SPECIES, BRACKETS );
    
    this->minimumDens2ResolvePPC = callPyFloatFunction( pInstance, GET_MIN_DENS_4_PPC, BRACKETS );
    
    this->maxTimestepsNum        = callPyFloatFunction( pInstance, GET_MAX_TIMESTEPS_NUM, BRACKETS );
    
    this->timestepsNum2Write     = callPyFloatFunction( pInstance, GET_TIMESTEP_WRITE, BRACKETS );
    
    this->outputDir              = callPyStringFunction( pInstance, GET_OUTPUT_DIR, BRACKETS );
    
    this->fileNameTemplate       = callPyStringFunction( pInstance, GET_FILENAME_TEMPLATE, BRACKETS );
    
    this->resistivity            = callPyFloatFunction( pInstance, GET_RESISTIVITY, BRACKETS );
    
    this->electronmass           = callPyFloatFunction( pInstance, GET_ELECTRON_MASS, BRACKETS );
    
    this->smoothStride     = (int) callPyLongFunction( pInstance, GET_PRESSURE_SMOOTH_STRIDE, BRACKETS);
    
    this->relaxFactor            = callPyFloatFunction( pInstance, GET_RELAX_FACTOR, BRACKETS );

    this->useIsothermalClosure = (int) callPyLongFunction( pInstance, IF2USE_ISOTHERMAL_CLOSURE, BRACKETS);
    this->electronTemperature  =       callPyFloatFunction( pInstance, GET_FIXED_ELECTRON_TEMPERATURE, BRACKETS );
    
    this->numOfSpots             = callPyFloatFunction( pInstance, NUMBER_OF_LASER_SPOTS, BRACKETS );
    
    
    if( numOfSpots > 0 ){
        
        numOfSpecies += 1; //need to reserve space for loaded particles
        
        this->prtclType2Load = callPyFloatFunction( pInstance, PARTICLE_TYPE2LOAD, BRACKETS );
        // -1 because in input file count starts in human manner from 1
        // prtclType2Load is needed to define mass and charge of the injected particles
        this->prtclType2Load = this->prtclType2Load - 1;
        
        if( prtclType2Load > (numOfSpecies-1) ){
            throw runtime_error("no such particle type for loading!");
        }
                
        this->laserPulseDuration_tsnum = (int) callPyLongFunction( pInstance, GET_LASER_PULSE_DURATION, BRACKETS);
    }
    
    auto callMethod = getPyMethod( pInstance, GET_CELL_BREAKDOWN_EFIELD_FACTOR, BRACKETS );
    
    if( callMethod != NULL ){
        this->cellBreakdownEfieldFactor = PyFloat_AsDouble(callMethod);
        string fctrMsg = "[Loader] [STABILITY] cell breakdown E-field factor = "
                                            +to_string(cellBreakdownEfieldFactor);
        logger.writeMsg(fctrMsg.c_str(), DEBUG);
    }else{
        this->cellBreakdownEfieldFactor = BIGN;
        string fctrMsg = "[Loader] [STABILITY] cellBreakdownEfieldFactor = "
                                    +to_string(cellBreakdownEfieldFactor)+" (default)";
        logger.writeMsg(fctrMsg.c_str(), DEBUG);
        
    }
    
    callMethod = getPyMethod( pInstance, GET_CRITICAL_PRESSURE, BRACKETS );
    
    if( callMethod != NULL ){
        this->criticalPressure = PyFloat_AsDouble(callMethod);
        string presMsg = "[Loader] [STABILITY] critical pressure = "
                                            +to_string(criticalPressure);
        logger.writeMsg(presMsg.c_str(), DEBUG);
    }else{
        this->criticalPressure = BIGN;
        string presMsg = "[Loader] [STABILITY] critical pressure = "
                            +to_string(criticalPressure)+" (default)";
        logger.writeMsg(presMsg.c_str(), DEBUG);
    }

    callMethod = getPyMethod( pInstance, GET_IONION_COLLISION_FREQUENCY_FACTOR, BRACKETS );
    if( callMethod != NULL ){
        this->collisionFrequencyFactor = PyFloat_AsDouble(callMethod);        
        string collMsg = "[Loader] [IONION_COLLISIONS] collisionFrequencyFactor = "
                                            +to_string(collisionFrequencyFactor);
        logger.writeMsg(collMsg.c_str(), DEBUG);
        callMethod = getPyMethod( pInstance, GET_DEFAULT_COLOUMB_LOG, BRACKETS );
        if( callMethod != NULL ){
            this->defaultCoulombLogarithm = PyFloat_AsDouble(callMethod);            
        }else{
            this->defaultCoulombLogarithm  = 25;
        }
        collMsg = "[Loader] [IONION_COLLISIONS] defaultCoulombLogarithm = "
                                            +to_string(defaultCoulombLogarithm);
        logger.writeMsg(collMsg.c_str(), DEBUG);
    }else{
        this->collisionFrequencyFactor = 0.0;
        this->defaultCoulombLogarithm  = 25;
    }
     
    
    
    
    if( rank == 0 ){

        string msg = runType == 0 ? "scratch" : "restart";
        msg = "[Loader] [COMMON] run from "+msg+" "+to_string(dim)+"D simulation";
        logger.writeMsg(msg.c_str(), INFO);
        
        msg = "[Loader] [COMMON]  Box size:";
        logger.writeMsg(msg.c_str(), INFO);

        for( int n = 0; n < 3; n++ ){
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
        msg = "[Loader] [COMMON] minimum density resolved by ppc number = "+to_string(minimumDens2ResolvePPC);
        logger.writeMsg(msg.c_str(), INFO);
        
        if( numOfSpots > 0 ){
            msg = "[Loader] [LASER] numOfSpots = "+to_string(numOfSpots)
                    +"; prtclType2Load = "+to_string(prtclType2Load+1);
            logger.writeMsg(msg.c_str(), INFO);
            msg = "[Loader] [LASER] laserPulseDuration_tsnum = "+to_string(laserPulseDuration_tsnum)
            +"; laserPulseDuration_omega = "+to_string(laserPulseDuration_tsnum*timeStep);
            logger.writeMsg(msg.c_str(), INFO);
        
        }
        
        msg = "[Loader] [OHM's LAW]: resistivity = "+to_string(resistivity);
        logger.writeMsg(msg.c_str(), INFO);

	    #ifdef USE_COLLISIONAL_RESIST_FACTOR
	       msg = "[Loader] [OHM's LAW]: temperature dependence (1/Te^1.5) factor for resistivity is ON";
           logger.writeMsg(msg.c_str(), INFO);
	    #endif
        
        #ifdef USE_EDGE_FACTOR
        msg = "[Loader] [OHM's LAW]: edge factor for ohm's law terms (pressure and hall) is ON";
        logger.writeMsg(msg.c_str(), INFO);
        #endif
        
	    if( useIsothermalClosure == 1){	
		  msg = "[Loader] [PRESSURE]: you use isothermal electrons, electron temperature is "+to_string(electronTemperature);
		  logger.writeMsg(msg.c_str(), INFO);
        }else{
        #ifdef IMPLICIT_PRESSURE
		msg = "[Loader] [PRESSURE]: you use implicit scheme ";
        #else
		msg = "[Loader] [PRESSURE]: you use explicit (subcycling) scheme (default)";
        #endif
		logger.writeMsg(msg.c_str(), INFO);
		msg = "[Loader] [PRESSURE]: relaxation factor = "+to_string(relaxFactor);
		logger.writeMsg(msg.c_str(), INFO);
        
		msg = "[Loader] [PRESSURE]: smoothing stride = "+to_string(smoothStride);
		logger.writeMsg(msg.c_str(), INFO);
        }
                
        logger.writeMsg("[Loader] initialization...OK!", DEBUG);
    }

}

double  Loader::getCollisionFrequencyFactor(){
    return collisionFrequencyFactor;
}

double  Loader::getDefaultCoulombLogarithm(){
    return defaultCoulombLogarithm;
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

vector<double> Loader::getVelocity( double x, double y, double z, int speciesType ){
    vector<double> velocity;
    
    for( int n = 0; n < 3; n++ ){
        string varName =GET_VELOCITY+dirs[n]+SPECIES+to_string(speciesType+1);

        double val = callPyFloatFunctionWith3args( pInstance, varName, BRACKETS_3DOUBLE, x, y, z );
        velocity.push_back(val);
    }
    return velocity;
}

vector<double> Loader::getVelocity4InjectedParticles( double x, double y, double z ){
    vector<double> velocity;
    
    for( int n = 0; n < 3; n++ ){
        string varName =GET_VELOCITY+dirs[n]+INJECTED_PARTICLES;

        double val = callPyFloatFunctionWith3args( pInstance, varName, BRACKETS_3DOUBLE, x, y, z );
        velocity.push_back(val);
    }
    return velocity;
}


vector<double> Loader::getFluidVelocity( double x, double y, double z, int speciesType ){
    vector<double> velocity;

    for( int n = 0; n < 3; n++ ){
        string varName =GET_FLUID_VELOCITY+dirs[n]+SPECIES+to_string(speciesType+1);

        double val = callPyFloatFunctionWith3args( pInstance, varName, BRACKETS_3DOUBLE, x, y, z );
        velocity.push_back(val);
    }
    return velocity;
}

vector<double> Loader:: getFluidVelocity4InjectedParticles( double x, double y, double z ){
    vector<double> velocity;

    for( int n = 0; n < 3; n++ ){
        string varName = GET_FLUID_VELOCITY+dirs[n]+INJECTED_PARTICLES;

        double val = callPyFloatFunctionWith3args( pInstance, varName, BRACKETS_3DOUBLE, x, y, z );
        velocity.push_back(val);
    }
    return velocity;
}

vector<double> Loader::getBfield( double x, double y, double z ){
    vector<double> bField;
    for( int n = 0; n < 3; n++ ){
        string varName = GET_BFIELD+dirs[n];
        double val = callPyFloatFunctionWith3args( pInstance, varName, BRACKETS_3DOUBLE, x, y, z );
        bField.push_back(val);
    }
    return bField;
}

double Loader::getDensity( double x, double y, double z, int speciesType ){
    string varName = GET_DENSITY+SPECIES+to_string(speciesType+1);
    return callPyFloatFunctionWith3args( pInstance, varName, BRACKETS_3DOUBLE, x, y, z );
}

double Loader::getMass4species( int speciesType ){
    string varName;
    if( numOfSpots > 0 && speciesType == numOfSpecies-1 ){
        varName = GET_MASS+SPECIES+to_string(prtclType2Load+1);
    }else{
        varName = GET_MASS+SPECIES+to_string(speciesType+1);
    }
    return callPyFloatFunction( pInstance, varName, BRACKETS );
}

double Loader::getCharge4species( int speciesType ){
    string varName;
    if( numOfSpots > 0 && speciesType == numOfSpecies-1 ){
        varName = GET_CHARGE+SPECIES+to_string(prtclType2Load+1);
    }else{
        varName = GET_CHARGE+SPECIES+to_string(speciesType+1);
    }
    return callPyFloatFunction( pInstance, varName, BRACKETS );
}

int Loader::getIfSpeciesFrozen( int speciesType ){
    if( numOfSpots > 0 && speciesType == numOfSpecies-1){
        return 0;// loaded particles are always unfrozen / needed for restart
    }
    string varName = GET_IFPARTICLETYPEISFROZEN+SPECIES+to_string(speciesType+1);
    return int(callPyFloatFunction( pInstance, varName, BRACKETS ));
}

int Loader::getDFtype(int speciesType){
    string varName = GET_DFTYPE+SPECIES+to_string(speciesType+1);
    return int(callPyFloatFunction( pInstance, varName, BRACKETS ));
}

int Loader::getDFtype4InjectedParticles(){
    string varName = GET_DFTYPE4INJECTEDPARTICLES;
    return int(callPyFloatFunction( pInstance, varName, BRACKETS ));
}


double Loader::getPPC4species(int speciesType){
    string varName;
    
    if( numOfSpots > 0 && speciesType == numOfSpecies-1 ){
         varName = GET_PPC4LOADED_PARTICLES;
    }else{
        varName = GET_PPC+SPECIES+to_string(speciesType+1);
    }
   
    return callPyFloatFunction( pInstance, varName, BRACKETS );
}


double Loader::getElectronPressureXX( double x, double y, double z ){
    string varName = GET_ELEPRESXX;
    return callPyFloatFunctionWith3args( pInstance, varName, BRACKETS_3DOUBLE, x, y, z );
}


double Loader::getElectronPressureYY( double x, double y, double z ){
    string varName = GET_ELEPRESYY;
    return callPyFloatFunctionWith3args( pInstance, varName, BRACKETS_3DOUBLE, x, y, z );
}


double Loader::getElectronPressureZZ( double x, double y, double z ){
    string varName = GET_ELEPRESZZ;
    return callPyFloatFunctionWith3args( pInstance, varName, BRACKETS_3DOUBLE, x, y, z );
}


double Loader::getTargetIonDensityProfile( double x, double y, double z ){
    string varName = GET_TARGET_ION_DENSITY2SUSTAIN;
    return callPyFloatFunctionWith3args( pInstance, varName, BRACKETS_3DOUBLE, x, y, z );
}

double Loader::getElectronPressureProfile( double x, double y, double z ){
    string varName = GET_ELECTRON_PRESSURE2SUSTAIN;
    return callPyFloatFunctionWith3args( pInstance, varName, BRACKETS_3DOUBLE, x, y, z );
}


Loader::~Loader(){
    Py_Finalize();
    logger.writeMsg("[Loader] FINALIZE...OK!", DEBUG);
}
