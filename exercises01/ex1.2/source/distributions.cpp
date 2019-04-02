#include "distributions.h"
#include <cmath>

RandomDistributions::RandomDistributions() :
    Random()
{
}
RandomDistributions::~RandomDistributions(){
}

double RandomDistributions::Exponential( double lambda){
    double x=Rannyu();
    return -std::log(1-x)/lambda;
}

double RandomDistributions::Cauchy(double x0, double gamma){
    double x=Rannyu();
    return gamma*std::tan(M_PI*(x-0.5))+x0;
}
