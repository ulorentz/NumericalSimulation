#ifndef __DISTRIBUTIONS_H_
#define __DISTRIBUTIONS_H_

#include "random.h"

class RandomDistributions : public Random {
    public:
        //Inherithed methods
        /*
         *
         * void SetRandom(int * , int, int);
         * void SaveSeed();
         * double Rannyu(void);
         * double Rannyu(double min, double max);
         * double Gauss(double mean, double sigma);
         */
        RandomDistributions();
        ~RandomDistributions();

        //Generates a variable distributed as a cauchy distributions:
        // p_cauchy=lambda*(1/((x-x0)^2+lambda^2)/pi
        //Arguments:
        //  x0:     center of cauchy distribution
        //  gamma:  parameter of cauchy
        //
        //Returns:
        //  x distributed as p_cauchy
        double Cauchy(double x0, double gamma); 
      
      
        //Generates a variable distributed as an exponential distribution:
        // p_exp(y)=lambda exp(-lambda y)
        //Arguments:
        //  lambda: parameter of exponential
        //
        //Returns:
        //  x distributed as p_exp
        double Exponential( double lambda); 

       
};

#endif //_DISTRIBUTIONS_H_
