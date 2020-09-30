#ifndef Misc_hpp
#define Misc_hpp

#include <stdio.h>
#include <cmath>
#include <vector>
#include <string>

#define  IDX(i, j, k, n0, n1, n2) ((k)+(n2)*((j)+(n1)*((i)+(n0)*(0))))
#define  IDX_2d(i, j, k, n0, n1, n2) ((j)+(n1)*(i)+(0)*(k)+(0)*(n2)+(n0)*(0))


#define BIGN 1.0E10
#define EPS8 1.0E-8
#define EPS4 1.0E-4
#define EPSILON           (1.0E-18)
#define EPSILON_VOLUME    (1.0E-10 )
#define EPSILON_DOUBLE    (1.0E-36)
#define EPSILON_COMPARATOR    (1.0E-18)
#define PI           3.14159265358979323846 

#define RNM  ((double)rand()/(double)RAND_MAX)

#define   PREDICTOR   0
#define   CORRECTOR   1 


bool areSame(double, double);

void crossProd(double[3], double[3], double[3]);

double edgeProfile(double dens);


#endif /* Misc_hpp */
