#include "ClosureManager.hpp"

using namespace std;
using namespace chrono;


#define FORWARD 0
#define BACKWARD 1
#define ALPHA 0.5


ClosureManager::ClosureManager(std::shared_ptr<Loader> ldr,
                               std::shared_ptr<GridManager> gridMnr):
                               loader(move(ldr)), gridMgr(move(gridMnr)){
    logger.reset(new Logger());
    initialize();
    logger->writeMsg("[ClosureManager] create...OK", DEBUG);
}

ClosureManager::~ClosureManager(){
    // need for driver calculation
    delete[] driverNext;
    delete[] pressuNext;
    delete[] electrnVel;
    delete[] pressuInit;
    delete[] bfieldPrev;
    delete[] pressureDampingCoeff;
    
}

void ClosureManager::initialize(){
    int xRes = loader->resolution[0],
        yRes = loader->resolution[1],
        zRes = loader->resolution[2];
    int xResG2 = xRes+2, yResG2 = yRes+2, zResG2 = zRes+2;
    
    int nG2= xResG2*yResG2*zResG2;
    
    driverNext = new double[nG2*6*sizeof(double)];
    pressuNext = new double[nG2*6*sizeof(double)];
    pressuInit = new double[nG2*6*sizeof(double)];
    electrnVel = new double[nG2*3*sizeof(double)];
    
    emass = loader->electronmass;
    
    initPressure();
    initPressureDampingCoeff();
    
    int h, ijkG1, ijkG2;
   
    VectorVar** pressure = gridMgr->getVectorVariableOnG2(PRESSURE);
    double pres;
    
    for( ijkG2 = 0; ijkG2 < nG2; ijkG2++ ){
        for( h = 0; h < 6; h++ ){
            pres = pressure[ijkG2]->getValue()[h];
            pressuNext[6*ijkG2+h] = pres;
            driverNext[6*ijkG2+h] = pres;
        }
    }

    int xResG1 = xRes+1, yResG1 = yRes+1, zResG1 = zRes+1;
    
    int nG1= xResG1*yResG1*zResG1;
    
    bfieldPrev = new double[nG1*3*sizeof(double)];
    
    VectorVar** bFieldKeep = gridMgr->getVectorVariableOnG1(MAGNETIC);
    for( ijkG1 = 0; ijkG1 < nG1; ijkG1++ ){
        for( h = 0; h < 3; h++ ){
            bfieldPrev[3*ijkG1+h] = bFieldKeep[ijkG1]->getValue()[h];
        }
    }
    
    calculatePressure(PREDICTOR, -1);
    calculatePressure(CORRECTOR, -1);
}

void ClosureManager::initPressureDampingCoeff(){

    //needed for open boundary conditions
    int xRes = loader->resolution[0],
        yRes = loader->resolution[1],
        zRes = loader->resolution[2];
    int xResG2 = xRes+2, yResG2 = yRes+2, zResG2 = zRes+2;
    
    double domainShiftX = loader->boxCoordinates[0][0];
    double domainShiftY = loader->boxCoordinates[1][0];
    double domainShiftZ = loader->boxCoordinates[2][0];
    
    double dx = loader->spatialSteps[0];
    double dy = loader->spatialSteps[1];
    double dz = loader->spatialSteps[2];
    
    double Lx = loader->boxSizes[0],
           Ly = loader->boxSizes[1],
           Lz = loader->boxSizes[2];
    
    int idx2use;
    int nG2 = xResG2*yResG2*zResG2;
    int i, j, k, idx;
    
    pressureDampingCoeff = new double[nG2*sizeof(double)];
    
    for( idx = 0; idx < nG2; idx++ ){
        pressureDampingCoeff[idx] = 1.0;
    }
    
    double leftXWidth = loader->dampingBoundaryWidth[0][0];
    double rigtXWidth = loader->dampingBoundaryWidth[0][1];
    
    double leftYWidth = loader->dampingBoundaryWidth[1][0];
    double rigtYWidth = loader->dampingBoundaryWidth[1][1];
    
    double leftZWidth = loader->dampingBoundaryWidth[2][0];
    double rigtZWidth = loader->dampingBoundaryWidth[2][1];
    
    double normLengthXleft, normLengthXrigt,
           normLengthYleft, normLengthYrigt,
           normLengthZleft, normLengthZrigt, minCoef;
    
    double x, y, z;
    
    for( i = 0; i < xResG2; i++ ){
        for( j = 0; j < yResG2; j++ ){
            for( k = 0; k < zResG2; k++ ){
                
                x = i*dx + domainShiftX;
                y = j*dy + domainShiftY;
                z = k*dz + domainShiftZ;
                
                idx2use = IDX(i,j,k, xResG2,yResG2,zResG2);
                normLengthXleft = 1; normLengthXrigt = 1;
                normLengthYleft = 1; normLengthYrigt = 1;
                normLengthZleft = 1; normLengthZrigt = 1;
                minCoef = 1;
                
                
                if( leftXWidth > 0.0 ){
                    normLengthXleft = x/leftXWidth;
                    if( normLengthXleft >= 0 && normLengthXleft < minCoef ){
                        minCoef = normLengthXleft;
                    }
                }
                
                if( rigtXWidth > 0.0 ){
                    normLengthXrigt = ( Lx + dx - x )/rigtXWidth;
                    if( normLengthXrigt >= 0 && normLengthXrigt < minCoef ){
                        minCoef = normLengthXrigt;
                    }
                }
                
                if( leftYWidth > 0.0 ){
                    normLengthYleft = y/leftYWidth;
                    if( normLengthYleft >= 0 && normLengthYleft < minCoef){
                        minCoef = normLengthYleft;
                    }
                }
                
                if( rigtYWidth > 0.0 ){
                    normLengthYrigt = ( Ly + dy - y )/rigtYWidth;
                    if( normLengthYrigt >= 0 && normLengthYrigt < minCoef ){
                        minCoef = normLengthYrigt;
                    }
                }
                
                if( leftZWidth > 0.0 ){
                    normLengthZleft = z/leftZWidth;
                    if( normLengthZleft >= 0 && normLengthZleft < minCoef ){
                        minCoef = normLengthZleft;
                    }
                }
                
                if( rigtZWidth > 0.0 ){
                    normLengthZrigt = ( Lz + dz - z )/rigtZWidth;
                    if( normLengthZrigt >= 0 && normLengthZrigt < minCoef ){
                        minCoef = normLengthZrigt;
                    }
                }
                pressureDampingCoeff[idx2use] = minCoef;
            }
        }
    }
}


void ClosureManager::initPressure(){
    
    double pxx, pyy, pzz;
    
    int idx, idxOnG2;
    
    double domainShiftX = loader->boxCoordinates[0][0];
    double domainShiftY = loader->boxCoordinates[1][0];
    double domainShiftZ = loader->boxCoordinates[2][0];
    
    double dx = loader->spatialSteps[0];
    double dy = loader->spatialSteps[1];
    double dz = loader->spatialSteps[2];
    
    int xRes = loader->resolution[0];
    int yRes = loader->resolution[1];
    int zRes = loader->resolution[2];
    int xResG2 = xRes+2, yResG2 = yRes+2, zResG2 = zRes+2;
    double G2shift = 0.5;
    double x,y,z;
    int i,j,k,h;
    for( i = 0; i < xRes; i++){
        for( j = 0; j < yRes; j++) {
            for( k = 0; k < zRes; k++){
                
                idxOnG2 = IDX(i+1, j+1, k+1, xResG2, yResG2, zResG2);
                
                x = (i + G2shift)*dx + domainShiftX;
                y = (j + G2shift)*dy + domainShiftY;
                z = (k + G2shift)*dz + domainShiftZ;
        
                pxx = loader->getElectronPressureXX(x,y,z);
                pyy = loader->getElectronPressureYY(x,y,z);
                pzz = loader->getElectronPressureZZ(x,y,z);
        
                gridMgr->setVectorVariableForNodeG2(idxOnG2, PRESSURE, 0, pxx);
                gridMgr->setVectorVariableForNodeG2(idxOnG2, PRESSURE, 3, pyy);
                gridMgr->setVectorVariableForNodeG2(idxOnG2, PRESSURE, 5, pzz);
                
                gridMgr->setVectorVariableForNodeG2(idxOnG2, PRESSURE_AUX, 0, pxx);
                gridMgr->setVectorVariableForNodeG2(idxOnG2, PRESSURE_AUX, 3, pyy);
                gridMgr->setVectorVariableForNodeG2(idxOnG2, PRESSURE_AUX, 5, pzz);
        
            }
        }
    }
    
    gridMgr->sendBoundary2Neighbor(PRESSURE);
    gridMgr->sendBoundary2Neighbor(PRESSURE_AUX);
    
    gridMgr->applyBC(PRESSURE);
    gridMgr->applyBC(PRESSURE_AUX);

    VectorVar** pdriverr = gridMgr->getVectorVariableOnG2(DRIVER);
    int nG2= xResG2*yResG2*zResG2;
    const double* dr;
    
    for( idxOnG2 = 0; idxOnG2 < nG2; idxOnG2++ ){
        for( h = 0; h < 6; h++ ){
            dr = pdriverr[idxOnG2]->getValue();
            gridMgr->setVectorVariableForNodeG2(idxOnG2, DRIVER_AUX, h, dr[h]);
            pressuInit[6*idxOnG2+h] = dr[h];
        }
    }
}


void ClosureManager::calculatePressure(int phase, int i_time){

        if( loader->useIsothermalClosure == 1 ){
		calculateIsothermalPressure();
		return;
	}
    
#ifdef IMPLICIT_PRESSURE
        implicitPressure(phase, i_time) ;
#else
        subCycledPressure(phase, i_time) ;
#endif

}

void ClosureManager::calculateIsothermalPressure(){

    auto start_time = high_resolution_clock::now();
    
    int xRes = loader->resolution[0];
    int yRes = loader->resolution[1];
    int zRes = loader->resolution[2];
    int xResG2 = xRes+2, yResG2 = yRes+2, zResG2 = zRes+2;
    int nG2 = xResG2*yResG2*zResG2;

    double pressure = 0.0;
    double Tele = loader->electronTemperature;
    double dens = 0.0;

    VectorVar** density  = gridMgr->getVectorVariableOnG2(DENSELEC);

    for(int ijkG2 = 0; ijkG2 < nG2; ijkG2++ ){
	dens = density[ijkG2]->getValue()[0];
	pressure = dens*Tele;
        gridMgr->setVectorVariableForNodeG2(ijkG2, PRESSURE_SMO, 0, pressure);
	gridMgr->setVectorVariableForNodeG2(ijkG2, PRESSURE_SMO, 3, pressure);
	gridMgr->setVectorVariableForNodeG2(ijkG2, PRESSURE_SMO, 5, pressure);
    }

    auto end_time = high_resolution_clock::now();
    string msgend ="[ClosureManager] calculateIsothermalPressure() total duration = "+
    to_string(duration_cast<milliseconds>(end_time - start_time).count())+" ms";
    logger->writeMsg(msgend.c_str(), DEBUG);
}


void ClosureManager::setIsotropization(double pres[6], double iTerm[6]) {
    
    double trP = (pres[0]+pres[3]+pres[5])/3;
    double buf = loader->relaxFactor;
    
    iTerm[0] = -buf*(pres[0] - trP);
    iTerm[1] = -buf*pres[1];
    iTerm[2] = -buf*pres[2];
    iTerm[3] = -buf*(pres[3] - trP);
    iTerm[4] = -buf*pres[4];
    iTerm[5] = -buf*(pres[5] - trP);
}


void ClosureManager::transformMatrix(double old[3][3], double new_[3][3], double transit[3][3], int way) {
    int k, l, m, n;
    
    for( k = 0; k < 3; k++ ){
        for( l = 0; l < 3; l++ ){
            new_[k][l] = 0.0;
            
            for( m = 0; m < 3; m++ ){
                for( n = 0; n < 3; n++ ){
                    if( way == BACKWARD ){
                        new_[k][l] += transit[m][k]*old[m][n]*transit[n][l];
                    }else if( way == FORWARD ){
                        new_[k][l] += transit[k][m]*old[m][n]*transit[l][n];
                    }else{
                        throw runtime_error("no such a way to transform");
                    }
                }
            }
        }
    }
}



void ClosureManager::ortho( double bField[3], double aw[3][3] ){
    
    double bmod = sqrt(bField[0]*bField[0]+bField[1]*bField[1]+bField[2]*bField[2]);
    
    if( bmod < EPS8 ){
        aw[0][0] = 1.0; aw[0][1] = 0.0; aw[0][2] = 0.0;
        aw[1][0] = 0.0; aw[1][1] = 1.0; aw[1][2] = 0.0;
        aw[2][0] = 0.0; aw[2][1] = 0.0; aw[2][2] = 1.0;
    }else{
        aw[0][0] = bField[0]/bmod;
        aw[0][1] = bField[1]/bmod;
        aw[0][2] = bField[2]/bmod;
        
        if( aw[0][0] == aw[0][1] && aw[0][1] == aw[0][2] ){
            aw[1][0] =  aw[0][1];
            aw[1][1] = -aw[0][0];
            aw[1][2] =  0.0;
        } else {
            aw[1][0] = (aw[0][2] - aw[0][1]);
            aw[1][1] = (aw[0][0] - aw[0][2]);
            aw[1][2] = (aw[0][1] - aw[0][0]);
        }
        
        bmod = sqrt(aw[1][0]*aw[1][0] + aw[1][1]*aw[1][1] + aw[1][2]*aw[1][2]);
        aw[1][0] /= bmod;
        aw[1][1] /= bmod;
        aw[1][2] /= bmod;
        
        aw[2][0] = aw[0][1]*aw[1][2]-aw[0][2]*aw[1][1];
        aw[2][1] = aw[0][2]*aw[1][0]-aw[0][0]*aw[1][2];
        aw[2][2] = aw[0][0]*aw[1][1]-aw[0][1]*aw[1][0];
    }
}


void ClosureManager::implicitPressure(int phase, int i_time) {
    
    auto start_time = high_resolution_clock::now();
    
    int xRes = loader->resolution[0];
    int yRes = loader->resolution[1];
    int zRes = loader->resolution[2];
    int xResG2 = xRes+2, yResG2 = yRes+2, zResG2 = zRes+2;
    int xResG1 = xRes+1, yResG1 = yRes+1, zResG1 = zRes+1;

    int nG2 = xResG2*yResG2*zResG2;
    int nG1 = xResG1*yResG1*zResG1;

    int ijkG2, ijkG1, h, i, j, k, l, m;
    
    double pPrev[6];
    double cTerm[6];
    double iTerm[6];
    
    double omega;
    
    double vecB[3];
    double vecBnext[3];
    double unitB[3];
    double modulusB;
    
    double *pPrevAll = new double[nG2*6*sizeof(double)];
    
    double ts = loader->getTimeStep();
    
    string msg ="[ClosureManager] start to calculate Pressure: imlpicit scheme ";
    logger->writeMsg(msg.c_str(), DEBUG);

    VectorVar** pressure    = gridMgr->getVectorVariableOnG2(PRESSURE);
    VectorVar** pressureaux = gridMgr->getVectorVariableOnG2(PRESSURE_AUX);
    VectorVar** pdriverr    = gridMgr->getVectorVariableOnG2(DRIVER);
    
    int magField2use;
    switch (phase){
        case PREDICTOR:
            magField2use = MAGNETIC_AUX;
            for( ijkG2 = 0; ijkG2 < nG2; ijkG2++ ){
                for( h = 0; h < 6; h++ ){
                    pPrevAll[ijkG2*6+h] = pressureaux[ijkG2]->getValue()[h];
                }
            }
            break;
        case CORRECTOR:
            magField2use = MAGNETIC;
            for( ijkG2 = 0; ijkG2 < nG2; ijkG2++ ){
                for( h = 0; h < 6; h++ ){
                    pPrevAll[ijkG2*6+h] = pressure[ijkG2]->getValue()[h];
                }
            }
            break;
        default :
            throw runtime_error("no phase");
    }
    
    const int* neighbourhood = gridMgr->getNeighbourhoodOnG1();
    VectorVar** bField = gridMgr->getVectorVariableOnG1(magField2use);
    
    int neighbour, idxNeigbor;
    
    const double* dr;
    double rhs;
    
    double P[3][3];
    double F[3][3];
    double pXYZ[3][3];
    double fXYZ[3][3];
    double uvw[3][3];
    
    double kappa, k2, buff1, buff2, trP;
    
    for( i = 1; i < xRes + 1; i++ ){
        for( j = 1; j < yRes + 1; j++ ){
            for( k = 1; k < zRes + 1; k++ ){
                
                ijkG2 = IDX(i, j, k, xResG2, yResG2, zResG2);
                
                for( h = 0; h < 6; h++ ){
                    pPrev[h] = pPrevAll[ijkG2*6+h];
                }
                
                setIsotropization(pPrev, iTerm);
                
                ijkG1 = IDX(i, j, k, xResG1, yResG1, zResG1);
                
                for( neighbour = 0; neighbour < 8; neighbour++ ){
                    idxNeigbor = neighbourhood[8*ijkG1+neighbour];
                    for( h = 0; h < 3; h++ ){
                        vecB[h] += 0.125*bfieldPrev[3*idxNeigbor+h];
                    }
                }
                
                modulusB = sqrt(vecB[0]*vecB[0]+vecB[1]*vecB[1]+vecB[2]*vecB[2]);
                
                if( modulusB > EPS8 ){
                    
                    omega = modulusB/emass;
                    for( h = 0; h < 3; h++ ) {
                        unitB[h] = vecB[h]/modulusB;
                    }
                
                    cTerm[0] = -( 2.0*(pPrev[1]*unitB[2]-pPrev[2]*unitB[1]) );
                    
                    cTerm[1] = -( pPrev[2]*unitB[0]-pPrev[0]*unitB[2]
                                 +pPrev[3]*unitB[2]-pPrev[4]*unitB[1] );
                    
                    cTerm[2] = -( pPrev[0]*unitB[1]-pPrev[1]*unitB[0]
                                 +pPrev[4]*unitB[2]-pPrev[5]*unitB[1] );
                    
                    cTerm[3] = -( 2.0*(pPrev[4]*unitB[0]-pPrev[1]*unitB[2]) );
                    
                    cTerm[4] = -( pPrev[1]*unitB[1]-pPrev[3]*unitB[0]
                                 +pPrev[5]*unitB[0]-pPrev[2]*unitB[2] );
                    
                    cTerm[5] = -( 2.0*(pPrev[2]*unitB[1]-pPrev[4]*unitB[0]) );

                }else{
                    omega = 0.0;
                    for( h = 0; h < 6; h++ ){
                        cTerm[h] = 0.0;
                    }
                }
                
                dr = pdriverr[ijkG2]->getValue();
                h = 0;
                for( l = 0; l < 3; l++ ){
                    for( m = l; m < 3; m++ ){
                        rhs = dr[h] + (1-ALPHA)*omega*cTerm[h] + iTerm[h];
                        fXYZ[l][m] = pPrev[h] + ts*rhs*pressureDampingCoeff[ijkG2];
                        fXYZ[m][l] = fXYZ[l][m];
                        h++;
                    }
                }
                
                for( neighbour = 0; neighbour < 8; neighbour++ ){
                    idxNeigbor = neighbourhood[8*ijkG1+neighbour];
                    for( h = 0; h < 3; h++ ){
                        vecB[h] += 0.125*bField[idxNeigbor]->getValue()[h];
                    }
                }
                
                modulusB = sqrt(vecB[0]*vecB[0]+vecB[1]*vecB[1]+vecB[2]*vecB[2]);
                
                if( modulusB > EPS8 ){
                    omega = modulusB/emass;
                    for ( h = 0; h < 3; h++ ) {
                        unitB[h] = vecB[h]/modulusB;
                    }
                }else{
                    omega = 0.0;
                    unitB[0] = 1.0;
                    unitB[1] = 0.0;
                    unitB[2] = 0.0;
                }
                
                ortho(unitB, uvw);
                
                transformMatrix(fXYZ, F, uvw, FORWARD);
                
                
                if( fXYZ[0][0]>BIGN || fXYZ[0][1]>BIGN || fXYZ[0][2]>BIGN ||
                   fXYZ[1][1]>BIGN || fXYZ[1][2]>BIGN || fXYZ[2][2]>BIGN ||
                   F[0][0]>BIGN || F[0][1]>BIGN || F[0][2]>BIGN ||
                   F[1][1]>BIGN || F[1][2]>BIGN || F[2][2]>BIGN){
                
                    string err ="[ClosureManager] bom bom 1 \n fXYZ[0][0] = "+to_string(fXYZ[0][0])
                    +"\n fXYZ[0][1] = "+to_string(fXYZ[0][1])
                    +"\n fXYZ[0][2] = "+to_string(fXYZ[0][2])
                    +"\n fXYZ[1][1] = "+to_string(fXYZ[1][1])
                    +"\n fXYZ[1][2] = "+to_string(fXYZ[1][2])
                    +"\n fXYZ[2][2] = "+to_string(fXYZ[2][2])
                    +"\n F[0][1] = "+to_string(F[0][1])
                    +"\n F[0][1] = "+to_string(F[0][1])
                    +"\n F[0][2] = "+to_string(F[0][2])
                    +"\n F[1][1] = "+to_string(F[1][1])
                    +"\n F[1][2] = "+to_string(F[1][2])
                    +"\n F[2][2] = "+to_string(F[2][2])
                    +"\n   dr[0] = "+to_string(dr[0])
                    +"\n   dr[1] = "+to_string(dr[1])
                    +"\n   dr[2] = "+to_string(dr[2])
                    +"\n   dr[3] = "+to_string(dr[3])
                    +"\n   dr[4] = "+to_string(dr[4])
                    +"\n   dr[5] = "+to_string(dr[5])
                    +"\n    cTerm[0] = "+to_string(cTerm[0])
                    +"\n    cTerm[1] = "+to_string(cTerm[1])
                    +"\n    cTerm[2] = "+to_string(cTerm[2])
                    +"\n    cTerm[3] = "+to_string(cTerm[3])
                    +"\n    cTerm[4] = "+to_string(cTerm[4])
                    +"\n    cTerm[5] = "+to_string(cTerm[5])
                    +"\n  Bunit[0] = "+to_string(unitB[0])
                    +"\n  Bunit[1] = "+to_string(unitB[1])
                    +"\n  Bunit[2] = "+to_string(unitB[2])
                    +"\n  uvw[0][0] = "+to_string(uvw[0][0])
                    +"\n  uvw[0][1] = "+to_string(uvw[0][1])
                    +"\n  uvw[0][2] = "+to_string(uvw[0][2])
                    +"\n  uvw[1][0] = "+to_string(uvw[1][0])
                    +"\n  uvw[1][1] = "+to_string(uvw[1][1])
                    +"\n  uvw[1][2] = "+to_string(uvw[1][2])
                    +"\n  uvw[2][0] = "+to_string(uvw[2][0])
                    +"\n  uvw[2][1] = "+to_string(uvw[2][1])
                    +"\n  uvw[2][2] = "+to_string(uvw[2][2])
                    +"\n omega = "+to_string(omega)
                    +"\n   iTerm[0] = "+to_string(iTerm[0])
                    +"\n   iTerm[1] = "+to_string(iTerm[1])
                    +"\n   iTerm[2] = "+to_string(iTerm[2])
                    +"\n   iTerm[3] = "+to_string(iTerm[3])
                    +"\n   iTerm[4] = "+to_string(iTerm[4])
                    +"\n   iTerm[5] = "+to_string(iTerm[5])
                    +"\n  ijkG2 = "+to_string(ijkG2)
                    +"\n i = "+to_string(i)
                    +"\n j = "+to_string(j)
                    +"\n k = "+to_string(k)
                    +"\n   +++++++++++++++++++++++++++++++++++";;
                    logger->writeMsg(err.c_str(), CRITICAL);
                    exit(-1);
                }
                
                kappa = ALPHA*omega*ts;
                
                k2 = kappa*kappa;
                buff1 = 1.0+k2;
                buff2 = 1.0+4.0*k2;
                
                P[0][0] = F[0][0];
                P[0][1] = (F[0][1]-kappa*F[0][2])/buff1;
                P[0][2] = (F[0][2]+kappa*F[0][1])/buff1;
                P[1][0] = P[0][1];
                P[2][0] = P[0][2];
                P[1][1] = (F[1][1]*(1.0+2.0*k2)-2.0*kappa*F[1][2]+2.0*k2*F[2][2])/buff2;
                P[1][2] = (kappa*F[1][1]+F[1][2]-kappa*F[2][2])/buff2;
                P[2][1] = P[1][2];
                P[2][2] = (2.0*k2*F[1][1]+2.0*kappa*F[1][2]+(1.0+2.0*k2)*F[2][2])/buff2;
                
                transformMatrix(P, pXYZ, uvw, BACKWARD);
                
                if( pXYZ[0][0]>BIGN || pXYZ[0][1]>BIGN || pXYZ[0][2]>BIGN ||
                   pXYZ[1][1]>BIGN || pXYZ[1][2]>BIGN || pXYZ[2][2]>BIGN ||
                   P[0][0]>BIGN || P[0][1]>BIGN || P[0][2]>BIGN ||
                   P[1][1]>BIGN || P[1][2]>BIGN || P[2][2]>BIGN){
                    
                    string err ="[ClosureManager] bom bom 2 \n pXYZ[0][0] = "+to_string(pXYZ[0][0])
                    +"\n pXYZ[0][1] = "+to_string(pXYZ[0][1])
                    +"\n pXYZ[0][2] = "+to_string(pXYZ[0][2])
                    +"\n pXYZ[1][1] = "+to_string(pXYZ[1][1])
                    +"\n pXYZ[1][2] = "+to_string(pXYZ[1][2])
                    +"\n pXYZ[2][2] = "+to_string(pXYZ[2][2])
                    +"\n P[0][1] = "+to_string(P[0][1])
                    +"\n P[0][1] = "+to_string(P[0][1])
                    +"\n P[0][2] = "+to_string(P[0][2])
                    +"\n P[1][1] = "+to_string(P[1][1])
                    +"\n P[1][2] = "+to_string(P[1][2])
                    +"\n P[2][2] = "+to_string(P[2][2])
                    +"\n   dr[0] = "+to_string(dr[0])
                    +"\n   dr[1] = "+to_string(dr[1])
                    +"\n   dr[2] = "+to_string(dr[2])
                    +"\n   dr[3] = "+to_string(dr[3])
                    +"\n   dr[4] = "+to_string(dr[4])
                    +"\n   dr[5] = "+to_string(dr[5])
                    +"\n    cTerm[0] = "+to_string(cTerm[0])
                    +"\n    cTerm[1] = "+to_string(cTerm[1])
                    +"\n    cTerm[2] = "+to_string(cTerm[2])
                    +"\n    cTerm[3] = "+to_string(cTerm[3])
                    +"\n    cTerm[4] = "+to_string(cTerm[4])
                    +"\n    cTerm[5] = "+to_string(cTerm[5])
                    +"\n  Bunit[0] = "+to_string(unitB[0])
                    +"\n  Bunit[1] = "+to_string(unitB[1])
                    +"\n  Bunit[2] = "+to_string(unitB[2])
                    +"\n  uvw[0][0] = "+to_string(uvw[0][0])
                    +"\n  uvw[0][1] = "+to_string(uvw[0][1])
                    +"\n  uvw[0][2] = "+to_string(uvw[0][2])
                    +"\n  uvw[1][0] = "+to_string(uvw[1][0])
                    +"\n  uvw[1][1] = "+to_string(uvw[1][1])
                    +"\n  uvw[1][2] = "+to_string(uvw[1][2])
                    +"\n  uvw[2][0] = "+to_string(uvw[2][0])
                    +"\n  uvw[2][1] = "+to_string(uvw[2][1])
                    +"\n  uvw[2][2] = "+to_string(uvw[2][2])
                    +"\n omega = "+to_string(omega)
                    +"\n   iTerm[0] = "+to_string(iTerm[0])
                    +"\n   iTerm[1] = "+to_string(iTerm[1])
                    +"\n   iTerm[2] = "+to_string(iTerm[2])
                    +"\n   iTerm[3] = "+to_string(iTerm[3])
                    +"\n   iTerm[4] = "+to_string(iTerm[4])
                    +"\n   iTerm[5] = "+to_string(iTerm[5])
                    +"\n  ijkG2 = "+to_string(ijkG2)
                    +"\n i = "+to_string(i)
                    +"\n j = "+to_string(j)
                    +"\n k = "+to_string(k)
                    +"\n   +++++++++++++++++++++++++++++++++++";;
                    logger->writeMsg(err.c_str(), CRITICAL);
                    exit(-1);
                }
                h = 0;
                for( l = 0; l < 3; l++ ){
                    for( m = l; m < 3; m++ ){
                        pPrevAll[ijkG2*6+h] = P[l][m];
                        h++;
                    }
                }
            }
        }
    }
    
    auto end_time1 = high_resolution_clock::now();
    string msg1 ="[ClosureManager] calculatePressure(): implicit duration = "+
                    to_string(duration_cast<milliseconds>(end_time1 - start_time).count())+" ms";
    logger->writeMsg(msg1.c_str(), DEBUG);
    
    auto end_time2 = high_resolution_clock::now();
    
    for( i = 1; i < xRes + 1; i++ ){
        for( j = 1; j < yRes + 1; j++ ){
            for( k = 1; k < zRes + 1; k++ ){
                ijkG2 = IDX(i, j, k, xResG2, yResG2, zResG2);
                for( h = 0; h < 6; h++ ){
                     gridMgr->setVectorVariableForNodeG2(ijkG2, PRESSURE, h, pPrevAll[ijkG2*6+h]);
                }
            }
        }
    }
    delete[] pPrevAll;
    
    gridMgr->sendBoundary2Neighbor(PRESSURE);
    gridMgr->applyBC(PRESSURE);
    
    if( i_time % loader->smoothStride == 0 ){
        gridMgr->smooth(PRESSURE);
        gridMgr->applyBC(PRESSURE);
    }
    
    for( ijkG2 = 0; ijkG2 < nG2; ijkG2++ ){
        for( h = 0; h < 6; h++) {
            gridMgr->setVectorVariableForNodeG2(ijkG2, PRESSURE_SMO, h,
                                                pressure[ijkG2]->getValue()[h]);
        }
    }
    
    gridMgr->smooth(PRESSURE_SMO);
    gridMgr->applyBC(PRESSURE_SMO);
    
    for( ijkG2 = 0; ijkG2 < nG2; ijkG2++ ){
        for( h = 0; h < 6; h++ ){
            pressuNext[6*ijkG2+h] = pressure[ijkG2]->getValue()[h];
        }
    }
    
    if( phase == PREDICTOR ) {
        for( ijkG2 = 0; ijkG2 < nG2; ijkG2++ ){
            for( h = 0; h < 6; h++ ){
                gridMgr->setVectorVariableForNodeG2(ijkG2, PRESSURE_AUX, h,
                                                    pressure[ijkG2]->getValue()[h]);
            }
        }
        
        VectorVar** bFieldKeep = gridMgr->getVectorVariableOnG1(MAGNETIC_AUX);
        for( ijkG1 = 0; ijkG1 < nG1; ijkG1++ ){
            for( h = 0; h < 3; h++ ){
                bfieldPrev[3*ijkG1+h] = bFieldKeep[ijkG1]->getValue()[h];
            }
        }
    }

    auto end_time3 = high_resolution_clock::now();
    string msg23 ="[ClosureManager] calculatePressure() set values duration = "+
    to_string(duration_cast<milliseconds>(end_time3 - end_time2).count())+" ms";
    logger->writeMsg(msg23.c_str(), DEBUG);
    
    setDriver(phase);
    
    
    auto end_time4 = high_resolution_clock::now();
    string msg22 ="[ClosureManager] calculatePressure() set driver duration = "+
    to_string(duration_cast<milliseconds>(end_time4 - end_time3).count())+" ms";
    logger->writeMsg(msg22.c_str(), DEBUG);
    
    string msgend ="[ClosureManager] calculatePressure() total duration = "+
    to_string(duration_cast<milliseconds>(end_time4 - start_time).count())+" ms";
    logger->writeMsg(msgend.c_str(), DEBUG);
    
}


void ClosureManager::subCycledPressure(int phase, int i_time) {
    
    auto start_time = high_resolution_clock::now();
    
    int xRes = loader->resolution[0];
    int yRes = loader->resolution[1];
    int zRes = loader->resolution[2];
    int xResG2 = xRes+2, yResG2 = yRes+2, zResG2 = zRes+2;
    int xResG1 = xRes+1, yResG1 = yRes+1, zResG1 = zRes+1;
    
    int nG2 = xResG2*yResG2*zResG2;
    int nG1 = xResG1*yResG1*zResG1;
    
    double dx = loader->spatialSteps[0];
    double dy = loader->spatialSteps[1];
    double dz = loader->spatialSteps[2];
    
    double domainShiftX = loader->boxCoordinates[0][0];
    double domainShiftY = loader->boxCoordinates[1][0];
    double domainShiftZ = loader->boxCoordinates[2][0];
    
    int ijkG2, ijkG1, h, i, j, k, m;
    
    double pSub[6];
    double cTerm[6];
    double iTerm[6];
    
    double omega;
    
    double vecB[3];
    double vecBnext[3];
    double unitB[3];
    double modulusB;
    
    double *vecBstartAll = new double[nG2*3*sizeof(double)];
    double *vecBstepAll  = new double[nG2*3*sizeof(double)];
    
    double *pSubAll      = new double[nG2*6*sizeof(double)];
    double *iTermAll     = new double[nG2*6*sizeof(double)];;
    
    double ts = loader->getTimeStep();
    const double subDt = ts*emass;
    
    const int numOfSubStep = (int)(1.0/emass);
    
    string msg ="[ClosureManager] start to calculate Pressure: subcycling subDt = "
    +to_string(subDt)+" numOfSubStep = "+to_string(numOfSubStep);
    logger->writeMsg(msg.c_str(), DEBUG);
    
    VectorVar** pressure    = gridMgr->getVectorVariableOnG2(PRESSURE);
    VectorVar** pressure_smo = gridMgr->getVectorVariableOnG2(PRESSURE_SMO);
    VectorVar** pressureaux = gridMgr->getVectorVariableOnG2(PRESSURE_AUX);
    VectorVar** pdriverr    = gridMgr->getVectorVariableOnG2(DRIVER);
    
    int magField2use;
    switch (phase){
        case PREDICTOR:
            magField2use = MAGNETIC_AUX;
            
            for( ijkG2 = 0; ijkG2 < nG2; ijkG2++ ){
                for( h = 0; h < 6; h++ ){
                    pSubAll[ijkG2*6+h] = pressureaux[ijkG2]->getValue()[h];
                }
            }
            break;
        case CORRECTOR:
            magField2use = MAGNETIC;
            
            for( ijkG2 = 0; ijkG2 < nG2; ijkG2++ ){
                for( h = 0; h < 6; h++ ){
                    pSubAll[ijkG2*6+h] = pressure[ijkG2]->getValue()[h];
                }
            }
            break;
        default :
            throw runtime_error("no phase");
    }
    
    const int* neighbourhood = gridMgr->getNeighbourhoodOnG1();
    VectorVar** bField = gridMgr->getVectorVariableOnG1(magField2use);
    
    int neighbour, idxNeigbor;
    
    //precalculations
    for( i = 1; i < xRes + 1; i++ ){
        for( j = 1; j < yRes + 1; j++ ){
            for( k = 1; k < zRes + 1; k++ ){
                
                for( h = 0; h < 3; h++ ){
                    vecBnext[h] = 0.0;
                    vecB[h]     = 0.0;
                }
                
                ijkG1 = IDX(i, j, k, xResG1, yResG1, zResG1);
                
                for( neighbour = 0; neighbour < 8; neighbour++ ){
                    idxNeigbor = neighbourhood[8*ijkG1+neighbour];
                    for( h = 0; h < 3; h++ ){
                        vecBnext[h] += 0.125*bField[idxNeigbor]->getValue()[h];
                        vecB[h]     += 0.125*bfieldPrev[3*idxNeigbor+h];
                    }
                }
                
                ijkG2 = IDX(i, j, k, xResG2, yResG2, zResG2);
                
                for( h = 0; h < 3; h++ ){
                    vecBstartAll[ijkG2*3+h] = vecB[h];
                    vecBstepAll [ijkG2*3+h] = (vecBnext[h]-vecB[h])/numOfSubStep;
                }
                
                for( h = 0; h < 6; h++ ){
                    pSub[h] = pSubAll[ijkG2*6+h];
                }
                
                setIsotropization(pSub, iTerm);
                
                for( h = 0; h < 6; h++ ){
                    iTermAll[ijkG2*6+h] = iTerm[h];
                }
                
            }
        }
    }
    
    auto end_time1 = high_resolution_clock::now();
    string msg1 ="[ClosureManager] calculatePressure(): precalculations duration = "+
    to_string(duration_cast<milliseconds>(end_time1 - start_time).count())+" ms";
    logger->writeMsg(msg1.c_str(), DEBUG);
    
    const double* dr;
    double trP, rhs;
    
    for( i = 1; i < xRes + 1; i++ ){
        for( j = 1; j < yRes + 1; j++ ){
            for( k = 1; k < zRes + 1; k++ ){
                
                ijkG2 = IDX(i, j, k, xResG2, yResG2, zResG2);
                
                for( h = 0; h < 6; h++ ){
                    pSub[h] = pSubAll[ijkG2*6+h];
                }
                
                dr = pdriverr[ijkG2]->getValue();
                
                for( m = 0; m < numOfSubStep; m++ ){
                    
                    for( h = 0; h < 6; h++ ){
                        cTerm[h] = 0.0;
                    }
                    
                    for( h = 0; h < 3; h++ ){
                        vecB[h] = vecBstartAll[ijkG2*3+h] + m*vecBstepAll[ijkG2*3+h];
                        unitB[h] = 0.0;
                    }
                    
                    modulusB = sqrt(vecB[0]*vecB[0]+vecB[1]*vecB[1]+vecB[2]*vecB[2]);
                    
                    omega = 0.0;
                    
                    if( modulusB > EPS8 ){
                        
                        for( h = 0; h < 3; h++ ){
                            unitB[h] = vecB[h]/modulusB;
                        }
                        omega = modulusB/emass;
                        
                        cTerm[0] = -( 2.0*(pSub[1]*unitB[2]-pSub[2]*unitB[1]) );
                        
                        cTerm[1] = -( pSub[2]*unitB[0]-pSub[0]*unitB[2]
                                     +pSub[3]*unitB[2]-pSub[4]*unitB[1] );
                        
                        cTerm[2] = -( pSub[0]*unitB[1]-pSub[1]*unitB[0]
                                     +pSub[4]*unitB[2]-pSub[5]*unitB[1] );
                        
                        cTerm[3] = -( 2.0*(pSub[4]*unitB[0]-pSub[1]*unitB[2]) );
                        
                        cTerm[4] = -( pSub[1]*unitB[1]-pSub[3]*unitB[0]
                                     +pSub[5]*unitB[0]-pSub[2]*unitB[2] );
                        
                        cTerm[5] = -( 2.0*(pSub[2]*unitB[1]-pSub[4]*unitB[0]) );
                    }
                    double damp = pressureDampingCoeff[ijkG2];
                    for ( h = 0; h < 6; h++ ){
                        rhs = dr[h] + omega*cTerm[h] + iTermAll[ijkG2*6+h];
                        pSub[h] += subDt*rhs*damp;
                    }
    
                }
                
                trP = pSub[0]+pSub[3]+pSub[5];
               if( trP < 0.0 || pSub[0] < 0 || pSub[3] < 0 || pSub[5] < 0 ){
                   
#ifdef HEAVYLOG
                   string err ="[ClosureManager] negative temperature: \n pSub[0] = "+to_string(pSub[0])
                   +"\n pSub[1] = "+to_string(pSub[1])
                   +"\n pSub[2] = "+to_string(pSub[2])
                   +"\n pSub[3] = "+to_string(pSub[3])
                   +"\n pSub[4] = "+to_string(pSub[4])
                   +"\n pSub[5] = "+to_string(pSub[5])
                   +"\n   dr[0] = "+to_string(dr[0])
                   +"\n   dr[1] = "+to_string(dr[1])
                   +"\n   dr[2] = "+to_string(dr[2])
                   +"\n   dr[3] = "+to_string(dr[3])
                   +"\n   dr[4] = "+to_string(dr[4])
                   +"\n   dr[5] = "+to_string(dr[5])
                   +"\n    cTerm[0] = "+to_string(cTerm[0])
                   +"\n    cTerm[1] = "+to_string(cTerm[1])
                   +"\n    cTerm[2] = "+to_string(cTerm[2])
                   +"\n    cTerm[3] = "+to_string(cTerm[3])
                   +"\n    cTerm[4] = "+to_string(cTerm[4])
                   +"\n    cTerm[5] = "+to_string(cTerm[5])
                   +"\n  Bunit[0] = "+to_string(unitB[0])
                   +"\n  Bunit[1] = "+to_string(unitB[1])
                   +"\n  Bunit[2] = "+to_string(unitB[2])
                   +"\n omega = "+to_string(omega)
                   +"\n   iTermAll[0] = "+to_string(iTermAll[ijkG2*6+0])
                   +"\n   iTermAll[1] = "+to_string(iTermAll[ijkG2*6+1])
                   +"\n   iTermAll[2] = "+to_string(iTermAll[ijkG2*6+2])
                   +"\n   iTermAll[3] = "+to_string(iTermAll[ijkG2*6+3])
                   +"\n   iTermAll[4] = "+to_string(iTermAll[ijkG2*6+4])
                   +"\n   iTermAll[5] = "+to_string(iTermAll[ijkG2*6+5])
                   +"\n  ijkG2 = "+to_string(ijkG2)
                   +"\n  subcycle m = "+to_string(m)
                   +"\n x = "+to_string(domainShiftX + dx*i)
                   +"\n y = "+to_string(domainShiftY + dy*j)
                   +"\n z = "+to_string(domainShiftZ + dz*k)
                   +"\n i = "+to_string(i)
                   +"\n j = "+to_string(j)
                   +"\n k = "+to_string(k)
                   +"\n   Xshift = "+to_string(domainShiftX)
                   +"\n   Yshift = "+to_string(domainShiftY)
                   +"\n   Zshift = "+to_string(domainShiftZ)
                   +"\n   +++++++++++++++++++++++++++++++++++";;
                   logger->writeMsg(err.c_str(), CRITICAL);
#endif
                   
                    for( h = 0; h < 6; h++ ){
                        pSub[h] = pressuInit[6*ijkG2+h];
                    }
                }
                
                for( h = 0; h < 6; h++ ){
                    pSubAll[ijkG2*6+h] = pSub[h];
                }
            
            }
        }
    }
    
    auto end_time2 = high_resolution_clock::now();
    string msg212 ="[ClosureManager] calculatePressure() subcycled duration = "+
    to_string(duration_cast<milliseconds>(end_time2 - end_time1).count())+" ms";
    logger->writeMsg(msg212.c_str(), DEBUG);
    
    
    for( i = 1; i < xRes + 1; i++ ){
        for( j = 1; j < yRes + 1; j++ ){
            for( k = 1; k < zRes + 1; k++ ){
                ijkG2 = IDX(i, j, k, xResG2, yResG2, zResG2);
                for (h = 0; h < 6; h++) {
                    gridMgr->setVectorVariableForNodeG2(ijkG2, PRESSURE, h, pSubAll[ijkG2*6+h]);
                }
            }
        }
    }
    
    delete[] vecBstartAll;
    delete[] vecBstepAll;
    delete[] pSubAll;
    delete[] iTermAll;
    
    gridMgr->sendBoundary2Neighbor(PRESSURE);
    gridMgr->applyBC(PRESSURE);
    
    for( ijkG2 = 0; ijkG2 < nG2; ijkG2++ ){
        for( h = 0; h < 6; h++) {
            gridMgr->setVectorVariableForNodeG2(ijkG2, PRESSURE_SMO, h,
                                                pressure[ijkG2]->getValue()[h]);
        }
    }
    
    gridMgr->smooth(PRESSURE_SMO);
    gridMgr->applyBC(PRESSURE_SMO);
    
    
    for( ijkG2 = 0; ijkG2 < nG2; ijkG2++ ){
        for( h = 0; h < 6; h++ ){
            pressuNext[6*ijkG2+h] = pressure[ijkG2]->getValue()[h];
        }
    }
    
    if( phase == PREDICTOR ){
        for( ijkG2 = 0; ijkG2 < nG2; ijkG2++ ){
            for( h = 0; h < 6; h++ ){
                gridMgr->setVectorVariableForNodeG2(ijkG2, PRESSURE_AUX, h,
                                                    pressure[ijkG2]->getValue()[h]);
            }
        }
        
        VectorVar** bFieldKeep = gridMgr->getVectorVariableOnG1(MAGNETIC_AUX);
        for( ijkG1 = 0; ijkG1 < nG1; ijkG1++ ){
            for( h = 0; h < 3; h++ ){
                bfieldPrev[3*ijkG1+h] = bFieldKeep[ijkG1]->getValue()[h];
            }
        }
    }
    
    auto end_time3 = high_resolution_clock::now();
    string msg23 ="[ClosureManager] calculatePressure() set values duration = "+
    to_string(duration_cast<milliseconds>(end_time3 - end_time2).count())+" ms";
    logger->writeMsg(msg23.c_str(), DEBUG);
    
    setDriver(phase);
    
    
    auto end_time4 = high_resolution_clock::now();
    string msg22 ="[ClosureManager] calculatePressure() set driver duration = "+
    to_string(duration_cast<milliseconds>(end_time4 - end_time3).count())+" ms";
    logger->writeMsg(msg22.c_str(), DEBUG);
    
    string msgend ="[ClosureManager] calculatePressure() total duration = "+
    to_string(duration_cast<milliseconds>(end_time4 - start_time).count())+" ms";
    logger->writeMsg(msgend.c_str(), DEBUG);
    
}




void ClosureManager::setDriver(int phase){
    
    int xRes = loader->resolution[0],
        yRes = loader->resolution[1],
        zRes = loader->resolution[2];
    int xResG2 = xRes+2, yResG2 = yRes+2, zResG2 = zRes+2;
    int nG2 = xResG2*yResG2*zResG2;
    
    int ijkG2;
    
    int i, j, k, l, m, h;
    double dTerms[3][3];
    
    double *pDrive = new double[nG2*6*sizeof(double)];
    
    VectorVar** current = gridMgr->getVectorVariableOnG2(CURRENT);
    
    for( ijkG2 = 0; ijkG2 < nG2; ijkG2++ ){
        for( h = 0; h < 3; h++ ){
            gridMgr->setVectorVariableForNodeG2(ijkG2, CURRENT_AUX, h,
                                                current[ijkG2]->getValue()[h]);
        }
    }
    
    gridMgr->smooth(CURRENT_AUX);
    gridMgr->applyBC(CURRENT_AUX);
    
    VectorVar** current_aux = gridMgr->getVectorVariableOnG2(CURRENT_AUX);
    
    VectorVar** density  = gridMgr->getVectorVariableOnG2(DENSELEC);
    VectorVar** velocity = gridMgr->getVectorVariableOnG2(VELOCION);
    const double* jcur;
    const double* edens;
    const double* vion;
    
    double domainShiftX = loader->boxCoordinates[0][0];
    double domainShiftY = loader->boxCoordinates[1][0];
    double domainShiftZ = loader->boxCoordinates[2][0];
    
    double dx = loader->spatialSteps[0];
    double dy = loader->spatialSteps[1];
    double dz = loader->spatialSteps[2];
    double dl[3] = {dx, dy, dz};
    
    for( i = 0; i < xRes + 2; i++ ){
        for (j = 0; j < yRes + 2; j++ ){
            for (k = 0; k < zRes + 2; k++ ){
                
                ijkG2 = IDX(i, j, k, xResG2, yResG2, zResG2);
                
                jcur = current_aux[ijkG2]->getValue();
                edens = density[ijkG2]->getValue();
                vion = velocity[ijkG2]->getValue();
                double revertdens = edens[0] < EPS8 ? 0.0 : edgeProfile(edens[0])/edens[0];
                for( l = 0; l < 3; l++ ){
                    electrnVel[3*ijkG2+l] = vion[l]-(jcur[l]*revertdens);
                    
                    gridMgr->setVectorVariableForNodeG2(ijkG2, VELOCELE, l,
                                                electrnVel[3*ijkG2+l]);
                }
            }
        }
    }
    
    gridMgr->applyBC(VELOCELE);
    
    double pe[3][3];
    double nabV[3][3];
    double nabP[3][3][3];
    
    gradientsVelocity();
    
    VectorVar** velDiagDer  = gridMgr->getVectorVariableOnG2(DRIVER_DIAG);
    VectorVar** velCrossDer  = gridMgr->getVectorVariableOnG2(DRIVER_CROSS);
    
    int n, s;
    double divV;
    
    float upwindIDX[3][3];
    float sign;
    int diffIDX[3][2];
    
    for( i = 1; i < xRes + 1; i++ ){
        for( j = 1; j < yRes + 1; j++ ){
            for( k = 1; k < zRes + 1; k++ ){
                
                ijkG2 = IDX(i, j, k, xResG2, yResG2, zResG2);
                
                /* __ indices to make difference in each dir {x,y,z} __ */
                diffIDX[0][0] = IDX(i+1, j  , k  , xResG2, yResG2, zResG2);
                diffIDX[0][1] = IDX(i-1, j  , k  , xResG2, yResG2, zResG2);
                diffIDX[1][0] = IDX(i  , j+1, k  , xResG2, yResG2, zResG2);
                diffIDX[1][1] = IDX(i  , j-1, k  , xResG2, yResG2, zResG2);
                diffIDX[2][0] = IDX(i  , j  , k+1, xResG2, yResG2, zResG2);
                diffIDX[2][1] = IDX(i  , j  , k-1, xResG2, yResG2, zResG2);
                
                /* __ upwind coefficients defined by flow direction __ */
                for( l = 0; l < 3; l++ ){//Vele{x,y,z}
                    sign = copysignf(1, electrnVel[3*ijkG2+l]);
                    upwindIDX[l][0] = +0.5*(1.0-sign)/dl[l]; //sign=-1 - counterflow direction use (Pi-1 - Pi)
                    upwindIDX[l][1] =           sign /dl[l]; //center
                    upwindIDX[l][2] = -0.5*(1.0+sign)/dl[l]; //sign=+1 - flow direction use (Pi - Pi+1)
                }
                
                h = 0;
                for( l = 0; l < 3; l++ ){
                    for( m = l; m < 3; m++ ){
                        
                        pe[l][m] = pressuNext[6*ijkG2+h];
                        pe[m][l] = pe[l][m];
                        
                        for( s = 0; s < 3; s++ ){//d{x,y,z}  nabla P upwind scheme
                            nabP[s][l][m] = upwindIDX[s][0]*pressuNext[6*diffIDX[s][0]+h]
                            +upwindIDX[s][1]*pe[l][m]
                            +upwindIDX[s][2]*pressuNext[6*diffIDX[s][1]+h];
                            
                            nabP[s][m][l] = nabP[s][l][m];// by Pij symmetry
                        }
                        h++;
                    }
                }
                
                nabV[0][0] = velDiagDer[ijkG2]->getValue()[0];
                nabV[1][1] = velDiagDer[ijkG2]->getValue()[1];
                nabV[2][2] = velDiagDer[ijkG2]->getValue()[2];
                nabV[0][1] = velCrossDer[ijkG2]->getValue()[0];
                nabV[0][2] = velCrossDer[ijkG2]->getValue()[1];
                nabV[1][0] = velCrossDer[ijkG2]->getValue()[2];
                nabV[1][2] = velCrossDer[ijkG2]->getValue()[3];
                nabV[2][0] = velCrossDer[ijkG2]->getValue()[4];
                nabV[2][1] = velCrossDer[ijkG2]->getValue()[5];
                
                divV = nabV[0][0]+nabV[1][1]+nabV[2][2];
                
                for( l = 0; l < 3; l++ ){
                    for( m = l; m < 3; m++ ){
                        
                        /* __ P nabla . V __ */
                        dTerms[l][m] = -pe[l][m]*divV;
                        
                        for( n = 0; n < 3; n++ ){
                            /* __ V . nabla P __ */
                            dTerms[l][m] -= electrnVel[3*ijkG2+n]*nabP[n][l][m];
                            /* __ P . nabla V __ */
                            dTerms[l][m] -= pe[l][n]*nabV[n][m];
                            /* __ P . nabla V (transposed) __ */
                            dTerms[l][m] -= pe[m][n]*nabV[n][l];
                        }
                        /* __ & set symmetrical terms __ */
                        dTerms[m][l] = dTerms[l][m];

                    }
                }
                
                h = 0;
                for( l = 0; l < 3; l++ ){
                    for( m = l; m < 3; m++ ){
                        pDrive[ijkG2*6+h] = dTerms[l][m];
                        h++;
                    }
                }
            }
        }
    }
    
    VectorVar** driveaux = gridMgr->getVectorVariableOnG2(DRIVER_AUX);
    
    switch (phase) {
        case PREDICTOR:
            for( ijkG2 = 0; ijkG2 < nG2; ijkG2++ ){
                for (h = 0; h < 6; h++) {
                    driverNext[6*ijkG2+h] = -driverNext[6*ijkG2+h]+2.0*pDrive[ijkG2*6+h];
                    gridMgr->setVectorVariableForNodeG2(ijkG2, DRIVER    , h, driverNext[ijkG2*6+h]);
                    gridMgr->setVectorVariableForNodeG2(ijkG2, DRIVER_AUX, h, pDrive[ijkG2*6+h]);
                }
            }
            gridMgr->sendBoundary2Neighbor(DRIVER_AUX);
            gridMgr->applyBC(DRIVER_AUX);
            
            break;
        case CORRECTOR:
            for( ijkG2 = 0; ijkG2 < nG2; ijkG2++ ){
                for( h = 0; h < 6; h++ ){
                    driverNext[6*ijkG2+h] = 0.5*(pDrive[ijkG2*6+h]+driveaux[ijkG2]->getValue()[h]);
                    gridMgr->setVectorVariableForNodeG2(ijkG2, DRIVER, h, driverNext[ijkG2*6+h]);
                }
            }
            break;
        default:
            throw runtime_error("no phase");
    }
    
    gridMgr->sendBoundary2Neighbor(DRIVER);
    gridMgr->applyBC(DRIVER);
    
    delete[] pDrive;
}

    void ClosureManager::gradientsVelocity() {
        
        int xRes = loader->resolution[0];
        int yRes = loader->resolution[1];
        int zRes = loader->resolution[2];
        int xResG2 = xRes+2, yResG2 = yRes+2, zResG2 = zRes+2;
        
        double dx = loader->spatialSteps[0];
        double dy = loader->spatialSteps[1];
        double dz = loader->spatialSteps[2];
        double dl[3] = {dx, dy, dz};
        
        int zeroOrderNeighb[6][3]  =
           {{-1, 0, 0}, {+1, 0, 0},
            { 0,-1, 0}, { 0,+1, 0},
            { 0, 0,-1}, { 0, 0,+1}};
        
        int firstOrderNeighb[24][3] =
           {{-1,-1, 0}, {+1,-1, 0},
            {-1,+1, 0}, {+1,+1, 0},
            {-1, 0,-1}, {+1, 0,-1},
            {-1, 0,+1}, {+1, 0,+1},
            
            {-1,-1, 0}, {-1,+1, 0},
            {+1,-1, 0}, {+1,+1, 0},
            { 0,-1,-1}, { 0,+1,-1},
            { 0,-1,+1}, { 0,+1,+1},
            
            {-1, 0,-1}, {-1, 0,+1},
            {+1, 0,-1}, {+1, 0,+1},
            { 0,-1,-1}, { 0,-1,+1},
            { 0,+1,-1}, { 0,+1,+1}};
        
        const double k1 = 0.250;// 1/4
        const double k2 = 0.0625;// 1/16
        
        int idxRight, idxLeft;
        double nabV[3][3];
        int ijkG2, i, j, k, l, m, pairNum;
        
        for( i = 1; i < xRes + 1; i++ ){
            for( j = 1; j < yRes + 1; j++ ){
                for( k = 1; k < zRes + 1; k++ ){
                    
                    ijkG2 = IDX(i, j, k, xResG2, yResG2, zResG2);
                    
                    /* __  nabla V centered in space __ */
                    for( l = 0; l < 3; l++ ){//d{x,y,z}
                        for( m = 0; m < 3; m++ ){//Vele{x,y,z}
                            
                            idxRight = IDX(i+zeroOrderNeighb[2*l+1][0],
                                           j+zeroOrderNeighb[2*l+1][1],
                                           k+zeroOrderNeighb[2*l+1][2],
                                           xResG2,yResG2,zResG2);
                            
                            idxLeft =  IDX(i+zeroOrderNeighb[2*l+0][0],
                                           j+zeroOrderNeighb[2*l+0][1],
                                           k+zeroOrderNeighb[2*l+0][2],
                                           xResG2,yResG2,zResG2);
                            
                            nabV[l][m] = k1*(electrnVel[3*idxRight+m] - electrnVel[3*idxLeft+m]);
                            
                            for( pairNum = 0; pairNum < 4; pairNum++ ){
                                idxRight = IDX(i+firstOrderNeighb[8*l+2*pairNum+1][0],
                                               j+firstOrderNeighb[8*l+2*pairNum+1][1],
                                               k+firstOrderNeighb[8*l+2*pairNum+1][2],
                                               xResG2,yResG2,zResG2);
                                
                                idxLeft  = IDX(i+firstOrderNeighb[8*l+2*pairNum+0][0],
                                               j+firstOrderNeighb[8*l+2*pairNum+0][1],
                                               k+firstOrderNeighb[8*l+2*pairNum+0][2],
                                               xResG2,yResG2,zResG2);
                                // have to keep second order because of hedp applications
                                // where too hot plasma oscillations lead to cross derivatives (dxVy, dyVx, ...) overestimation
                                // and negative temperature as result
                                nabV[l][m] += k2*(electrnVel[3*idxRight+m] - electrnVel[3*idxLeft+m]);
                            }
                            nabV[l][m] = nabV[l][m]/dl[l];
                        }
                    }
                    
                    
                    gridMgr->setVectorVariableForNodeG2(ijkG2, DRIVER_DIAG, 0, nabV[0][0]);
                    gridMgr->setVectorVariableForNodeG2(ijkG2, DRIVER_DIAG, 1, nabV[1][1]);
                    gridMgr->setVectorVariableForNodeG2(ijkG2, DRIVER_DIAG, 2, nabV[2][2]);
                    
                    gridMgr->setVectorVariableForNodeG2(ijkG2, DRIVER_CROSS, 0, nabV[0][1]);
                    gridMgr->setVectorVariableForNodeG2(ijkG2, DRIVER_CROSS, 1, nabV[0][2]);
                    gridMgr->setVectorVariableForNodeG2(ijkG2, DRIVER_CROSS, 2, nabV[1][0]);
                    gridMgr->setVectorVariableForNodeG2(ijkG2, DRIVER_CROSS, 3, nabV[1][2]);
                    gridMgr->setVectorVariableForNodeG2(ijkG2, DRIVER_CROSS, 4, nabV[2][0]);
                    gridMgr->setVectorVariableForNodeG2(ijkG2, DRIVER_CROSS, 5, nabV[2][1]);
                }
            }
        }
        
        gridMgr->sendBoundary2Neighbor(DRIVER_DIAG);
        gridMgr->applyBC(DRIVER_DIAG);
        
        gridMgr->smooth(DRIVER_DIAG);
        gridMgr->applyBC(DRIVER_DIAG);
        
        gridMgr->sendBoundary2Neighbor(DRIVER_CROSS);
        gridMgr->applyBC(DRIVER_CROSS);
        
        gridMgr->smooth(DRIVER_CROSS);
        gridMgr->applyBC(DRIVER_CROSS);
        
    }
