#include "variationalMC.h"
#include <algorithm>        //std::min

Particle1D::Particle1D(double mu, 
                       double sigma, 
                       double start_x, 
                       double jump) :
    mu(mu),
    sigma(sigma),
    x(start_x),
    jump(jump),
    rand("Primes", "seed.in")
{
    in=0;
    tot=0;
    best_mu=NAN;
    best_sigma=NAN;
}


void Particle1D::ResetMetropolis(double mu, 
                                 double sigma,
                                 double start_x, 
                                 double jump)
{
    Particle1D::mu=mu;
    Particle1D::sigma=sigma;
    x=start_x;
    Particle1D::jump=jump;
    in=0;
    tot=0;
}

void Particle1D::FindOptimum(double starting_mu, 
                             double var_mu,
                             double starting_sigma,
                             double var_sigma,
                             double start_x,
                             double jump, 
                             unsigned int n_samples,
                             unsigned int n_integral, 
                             unsigned int iterations){

    double min_mu=0., min_sigma=0., sigma_try, mu_try,
           integral=1E10, integral_old=0., delta_mu, 
           delta_sigma;

    min_mu=starting_mu;
    min_sigma=starting_sigma;

    double tmp_mu;
    for(unsigned int i=0; i<iterations; i++){ //number of "sweep"
        delta_mu=2.*var_mu/(n_samples*std::pow(10.,i)); //first time is
                                                        //divided by nsamples
        delta_sigma=2.*var_sigma/(n_samples*std::pow(10.,i)); //same
        mu_try=min_mu-var_mu/std::pow(10.,i); //fist time it's just var_mu
        tmp_mu=mu_try;
        sigma_try=min_sigma-var_sigma/std::pow(10.,i); //same
        for(unsigned int k=0; k<n_samples; ++k){
            for(unsigned int z=0; z<n_samples; ++z){
                ResetMetropolis(mu_try, sigma_try, start_x, jump);
                for(unsigned int j=0; j<n_integral; ++j)
                    integral+=Sample();
    
                integral/=n_integral;
                if(integral < integral_old){
                    min_mu=mu_try; 
                    min_sigma=sigma_try;
                }
                integral_old=integral;
                integral=0.;
                mu_try+=delta_mu;
            }
            mu_try=tmp_mu;
            sigma_try+=delta_sigma;
        }
    }
    best_mu=min_mu;
    best_sigma=min_sigma;

}

void Particle1D::Metropolis(){
    x_try=(rand.Rannyu()-0.5)*jump*2.+x;
    if(rand.Rannyu()<std::min(1., Density(x_try)/Density(x))){
        in++;
        x=x_try;
    }
    tot++;
}


double Particle1D::Sample() {
    Metropolis();
    return x;
//    return F(x);
}
