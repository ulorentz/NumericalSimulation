#ifndef _VARIATIONAL_MC_H_
#define _VARIATIONAL_MC_H_
#include <cmath>
#include "random.h"

class Particle1D{
    public:
        //Class that is aimed to find the ground state of the following problem:
        //H(x)=p^2/2m + x^4+5*x^2/2
        //with h_bar and m both equal to unit.
        //And a trial wavefunction 
        //psi(x)_{sigma,mu}= N [exp{-(x-mu)^2/2sigma^2 +exp{-(x+mu)^2/2sigma^2}]
        //
        //The parameters of the constructor are:
        // mu: variational parameter of the wavefunction
        // sigma: variational parameter of the wavefunction
        // start_x: starting point for metropolis algorithm
        // jump: jump step for the metropolis algorithm
        Particle1D(double mu, 
                   double sigma, 
                   double start_x=0., 
                   double jump=1.);

        //Resets the metropolis algorithm. Parameters are the same defined in 
        //the class constructor.
        void ResetMetropolis(double mu, 
                             double sigma,
                             double start_x=0., 
                             double jump=1.);
        
        //samples and returns one point
        double Sample();
       

        double Energy(double x) const;

        //Returns the actual acceptance rate
        double Acceptance() const;


        // starting_mu: variational parameter of the wavefunction. The parameter 
        //              will be variated in the interval ->
        //              
        //              [starting_mu-var_mu, starting_mu+var_mu] 
        //              
        //              where var_mu is the following parameter:
        // var_mu: semi amplitude of the inteval where to sample mu.
        // starting_sigma: variational parameter of the wavefunction. 
        //                 The parameter will be variated in the interval ->
        //              
        //                 [starting_sigma-var_sigma, starting_sigma+var_sigma] 
        //              
        //                 where var_sigma is the following parameter:
        // var_signa: semi amplitude of the inteval where to sample sigma. 
        // start_x: starting point for metropolis algorithm
        // jump: jump step for the metropolis algorithm
        // n_samples: number of samples to try in the two above intervals
        // n_integral: the integral is computing as sum p(x_i) f(x_i), with
        //             i from one to N, n_integral is N.
        // iterations: number of interations for finding the optimum. At each 
        //             iteration the mu and sigma are centered in the optimum
        //             of previous iter, and vars are multiplied by 1/10.
        void FindOptimum(double starting_mu, 
                         double var_mu,
                         double starting_sigma,
                         double var_sigma,
                         double start_x,
                         double jump, 
                         unsigned int n_samples,
                         unsigned int n_integral, 
                         unsigned int iterations=1);

        //Used to get the best parameters found in "FindOptimum", pass as
        //reference the two parameters, self explaining names. 
        //If "FindOptimum" hasn't been run previously, both params will be set 
        //to "nan". 
        void GetOptimalParameters(double &mu, double &sigma);
    private:
        
        //density function defined as:
        // p(x)=|psi|^2/(integral{|psi|^2}) 
        double Density(double x) const;

        //metropolis move
        void Metropolis();

         //wavefunction parameters
        double mu, sigma;
        //actual metropolis point and trial one
        double x, x_try;
        //metropolis jump step
        double jump;
        //used for acceptance rate
        unsigned int in, tot;
    
        //optimal parameters
        double best_mu, best_sigma;    
        Random rand;
};


        
        

/******************/
/* Inline methods */
/******************/
//macro definitions because I'm lazy and because I don't want to waste memory 
//defining variables in the following functions
#define SIGMA (sigma*sigma)
#define POT(x) (std::pow(x,4)-5.*std::pow(x,2)/2.)
#define EXP1(x) std::exp(-(x-mu)*(x-mu)/(SIGMA*2.))
#define EXP2(x) std::exp(-(x+mu)*(x+mu)/(SIGMA*2.))
inline double Particle1D::Density(double x) const{
    return std::pow(EXP1(x)+EXP2(x),2);
}

inline double Particle1D::Energy(double x) const{
    double k=(x*x+mu*mu)/SIGMA-1.-(2.*x*mu/SIGMA)*std::tanh(x*mu/SIGMA);
    k/=-2.*SIGMA;
  return k+POT(x);
}

#undef SIGMA
#undef POT
#undef EXP1
#undef EXP2

inline double Particle1D::Acceptance() const{
    return double(in)/tot;
}

        
inline void Particle1D::GetOptimalParameters(double &mu, double &sigma){
    mu=best_mu;
    sigma=best_sigma;
}
#endif //_VARIATIONAL_MC_H_
