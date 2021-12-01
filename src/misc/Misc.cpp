
#include "Misc.hpp"
using namespace std;

bool areSame(double a, double b)
{
    return fabs(a - b) < EPSILON_COMPARATOR;
}



void crossProd(double A[3], double B[3], double cross[3]){
    cross[0] = A[1] * B[2] - A[2] * B[1];
    cross[1] = A[2] * B[0] - A[0] * B[2];
    cross[2] = A[0] * B[1] - A[1] * B[0];
}


double edgeProfile(double dens){
    
    double const densLow  = 0.000001;
    double const densHigh = 0.1;
    
    double x = (dens-densLow)/(densHigh-densLow)-1;
    
    if (x < -1.0) { x = -1.0; }
    else if (x > 0.0) { x = 0.0; }
    
    double modX = fabs(x);
    
    double res = -6.0*modX*modX*modX*modX*modX
    +15.0*x*x*x*x
    -10.0*modX*modX*modX
    +1;
    
#ifdef USE_EDGE_FACTOR
    return res;
#else
    return 1.0;
#endif
    
}

double edgeProfilePressure(double pres){
    
    double const presLow  = 0.001;
    double const presHigh = 10.0;
    
    double x = (pres-presLow)/(presHigh-presLow)-1;
    
    if (x < -1.0) { x = -1.0; }
    else if (x > 0.0) { x = 0.0; }
    
    double modX = fabs(x);
    
    double res = -6.0*modX*modX*modX*modX*modX
    +15.0*x*x*x*x
    -10.0*modX*modX*modX
    +1;
    
    return res;    
}


double polynomByRoch(double val){
    
    if (val < -1.0) { val = -1.0; }
    else if (val > 0.0) { val = 0.0; }
    
    double modX = fabs(val);
    
    double res = -6.0*modX*modX*modX*modX*modX
    +15.0*modX*modX*modX*modX
    -10.0*modX*modX*modX
    +1;
    
    return res;
}

