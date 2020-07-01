
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


