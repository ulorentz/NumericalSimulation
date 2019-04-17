#ifndef _HYDROGEN_H_
#define _HYDROGEN_H_
#include <vector>
#include <string>
#include "random.h"
#include "functions.h"

class hydrogen{
    public:
        //Arguments:
        //  - x y z are the starting point for the metropolis algorithm
        //  - step: is the step size in metropolis algorithm
        //  - probability_density: density to sample, implemented in this class 
        //    are the 100 hydrogenoid level, and the 210; available options, so,
        //    are "100" and "210".
        //  - random_backend: is the method for generating random number, 
        //    options are "uniform" and "gauss"
        hydrogen(double xstart, 
                 double ystart, 
                 double zstart, 
                 double step_size,
                 std::string probability_density="100",
                 std::string random_backend="uniform");
        
        ~hydrogen();
        
        //used to reset the starting point for the metropolis algorithm
        //Arguments:
        //  - x y z are the starting point for the metropolis algorithm
        //  - step: is the step size in metropolis algorithm
        //  - probability_density: density to sample, implemented in this class 
        //    are the 100 hydrogenoid level, and the 210; available options, so,
        //    are "100" and "210".
        //  - random_backend: is the method for generating random number, 
        //    options are "uniform" and "gauss"
        void reset_metropolis(double xstart, 
                              double ystart, 
                              double zstart, 
                              double step,
                              std::string probability_density="100",
                              std::string random_backend="uniform");

        //sample a new point with the probability density and backend random 
        //generator previously set. 
        void sample();

        //returns the actual x,y,z position previously sampled
        std::vector<double> get_position() const;
        
        //returns the actual radius of the position previously sampled
        double get_radius() const;

        //return fraction of accepted moves
        float acceptance() const;

    private:
        //private methods
        //metropolis methods
        //uniform
        void metropolis_uniform(); 
        //gauss
        void metropolis_gauss(); 
        
        //probability density for level 100
        double prob100(double x, double y, double z) const;
        //probability density for level 210
        double prob210(double x, double y, double z) const;
         
        //callable pointer with correct "backend" and prob density settings
        void (hydrogen::*run)();
        //callable pointer with set probability function
        double (hydrogen::*probability)(double x, double y, double z) const;

        
        //private variables
        double x, y, z;
        //for metropolis algorithm
        double xtry, ytry, ztry, alpha;
        double step;
        Random rand;
        unsigned int in, tot;
};  


// Inline Methods

inline float hydrogen::acceptance() const
{
    return float(in)/tot;
}

inline void hydrogen::sample()
{
    //calling the sampling method to wich run points to. 
    //Has been set in the constructor or in the reset method.
    (*this.*run)();
}

#endif //_HYDROGEN_H_
