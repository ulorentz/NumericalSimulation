#ifndef _ACCEPT_REJECT_H_
#define _ACCEPT_REJECT_H_

#include "random.h"
#include <cmath>


class AcceptReject : public Random {
public:
    //Inherits from Random the following methods. The seed must be initialized.
    //  void SetRandom(int * , int, int);
    //  void SaveSeed();
    //  double Rannyu(void);
    //  double Rannyu(double min, double max);
    //  double Gauss(double mean, double sigma);

    AcceptReject();
    ~AcceptReject();

    //Generates random numbers ditributed within [0, 1) as a specific function 
    //with the accept-reject method.
    //
    //Arguments:
    //  func:   is a pointer to a function that specifies the analytical form of 
    //          distribution. Note that the function MUST be positive defined 
    //          in the interval [0,1). For instance func could be:
    //              double func(double x){
    //                  return x*x;
    //              }
    //
    //  max:    spcifies the maximum value that function takes on [0,1), defualt
    //          is 1.
    //Returns:
    //  random number distributed as function in the interval [0,1).
    //
    //Example:
    //  with fun definited as above, if rand is an instance of AcceptReject, 
    //  should be called as:
    //      double x = rand.distribution(func, 1.);
    double distribution(double (*function)(double), double max=1);

};

#endif //_ACCEPT_REJECT_H_
