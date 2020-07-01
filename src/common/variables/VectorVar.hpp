#ifndef VectorVar_hpp
#define VectorVar_hpp

#include <stdio.h>
#include <vector>
#include <string>
#include <iostream>

class VectorVar{

private:
    int name;
    double value[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    int size = 6;
    
public:
    VectorVar();
    ~VectorVar();
    VectorVar(int, std::vector<double>);
    int getName();
    int getSize();
    const double* getValue();
    void setValue(std::vector<double>);
    void setValue(double*);
    void setValue(const double*);
    void setValue(int, double);
    void addValue(int, double);
};
#endif /* VectorVar_hpp */

