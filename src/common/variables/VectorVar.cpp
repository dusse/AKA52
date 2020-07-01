#include "VectorVar.hpp"

using namespace std;

VectorVar::VectorVar(){
}

VectorVar::VectorVar(int name, std::vector<double> valueNew){
    this->name = name;
    this->size = valueNew.size();
    for (int dim=0; dim<size; dim++){
        value[dim] = valueNew[dim];
    }
}


int VectorVar::getName(){
    return name;
}

const double* VectorVar::getValue(){
    return value;
}

int VectorVar::getSize(){
    return size;
}

void VectorVar::setValue(std::vector<double> newValue){
    for (int dim=0; dim<size; dim++){
        value[dim] = newValue[dim];
    }
}

void VectorVar::setValue(double* newValue){
    for (int dim=0; dim<size; dim++){
        value[dim] = newValue[dim];
    }
}

void VectorVar::setValue(const double* newValue){
    for (int dim=0; dim<size; dim++){
        value[dim] = newValue[dim];
    }
}

void VectorVar::setValue(int dim, double newValue){
        value[dim] = newValue;
}

void VectorVar::addValue(int dim, double newValue){
    value[dim] += newValue;
}

VectorVar::~VectorVar(){
}