#ifndef _DISCRETE_RANDOM_H_
#define _DISCRETE_RANDOM_H_
#include "random.h"

class DiscreteRandom : public Random{
    public:
        //Inherits from Random the following methods. The seed must be 
        //initialized.
        //  void SetRandom(int * , int, int);
        //  void SaveSeed();
        //  double Rannyu(void);
        //  double Rannyu(double min, double max);
        //  double Gauss(double mean, double sigma);

        DiscreteRandom();
        ~DiscreteRandom();

        //Method that generate a uniform natural (discrete) random variable in 
        //the interval [a, b] with the bounds included.
        //Arguments:
        //  int a: lower bound
        //  int b: upper bound
        //
        //Returns:
        //  int in the range [a, b].
        int uniform(int a, int b);
};

#endif //_DISCRETE_RANDOM_H_
