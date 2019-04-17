#include "hydrogen.h"
#include <algorithm>        //std::min
#include <cmath>
#include <iostream>

hydrogen::hydrogen(double xstart, 
                   double ystart, 
                   double zstart, 
                   double step_size,
                   std::string probability_density,
                   std::string random_backend):
    x(xstart),
    y(ystart),
    z(zstart),
    step(step_size)
{
    //not valid random backend
    if (random_backend!="uniform" and random_backend!="gauss"){
        std::cout << "Backend '" << random_backend << "' not found.\n"
            "Possibilities are: 'uniform', 'gauss'." << std::endl;
        exit(0);
    }
    //not valid probability density
    if (probability_density!="100" and probability_density!="210"){
        std::cout << "Probability density '"<< probability_density << "' not"
            "found.\nPossibilities are: '100', '210'.\n";
        exit(0);
    }
    if (random_backend=="uniform") {
        run=&hydrogen::metropolis_uniform;
        if (probability_density=="100")
            probability=&hydrogen::prob100;
        else
            probability=&hydrogen::prob210;
    }
    else{
        run=&hydrogen::metropolis_gauss;
        if (probability_density=="100")
            probability=&hydrogen::prob100;
        else
            probability=&hydrogen::prob210;
    }

    initialize_seed(rand, "Primes", "seed.in");

    in=0;
    tot=0;
}

hydrogen::~hydrogen()
{
}

void hydrogen::reset_metropolis(double xstart, 
                                double ystart, 
                                double zstart, 
                                double step_size,
                                std::string probability_density,
                                std::string random_backend)
{
    x=xstart;
    y=ystart;
    z=zstart;
    step=step_size;
    //not valid random backend
    if (random_backend!="uniform" and random_backend!="gauss"){
        std::cout << "Backend '" << random_backend << "' not found.\n"
            "Possibilities are: 'uniform', 'gauss'." << std::endl;
        exit(0);
    }
    //not valid probability density
    if (probability_density!="100" and probability_density!="210"){
        std::cout << "Probability density '"<< probability_density << "' not"
            "found.\nPossibilities are: '100', '210'.\n";
        exit(0);
    }
    if (random_backend=="uniform") {
        run=&hydrogen::metropolis_uniform;
        if (probability_density=="100")
            probability=&hydrogen::prob100;
        else
            probability=&hydrogen::prob210;
    }
    else{
        run=&hydrogen::metropolis_gauss;
        if (probability_density=="100")
            probability=&hydrogen::prob100;
        else
            probability=&hydrogen::prob210;
    }

    in=0;
    tot=0;
}

double hydrogen::prob100(double x_, double y_, double z_) const
{
    //working on bohr units
    return std::exp(-2.*std::sqrt(x_*x_+y_*y_+z_*z_))/M_PI;
}

double hydrogen::prob210(double x_, double y_, double z_) const
{
    double rsquare=x_*x_+y_*y_+z_*z_;
    //working on bohr units
    return std::exp(-std::sqrt(rsquare))*x_*x_/(32.*M_PI);
}


void hydrogen::metropolis_uniform()
{
    xtry=x+(rand.Rannyu()-0.5)*step; 
    ytry=y+(rand.Rannyu()-0.5)*step; 
    ztry=z+(rand.Rannyu()-0.5)*step; 
    alpha=std::min(1., (*this.*probability)(xtry, ytry, ztry)/
            (*this.*probability)(x,y,z));
    if (rand.Rannyu()<=alpha) { //valid even if alpha = 1
        x=xtry;
        y=ytry;
        z=ztry;
        in++;
    }
    tot++;
}

void hydrogen::metropolis_gauss()
{
    xtry=rand.Gauss(x,step/2.);
    ytry=rand.Gauss(y,step/2.);
    ztry=rand.Gauss(z,step/2.);
    alpha=std::min(1., (*this.*probability)(xtry, ytry, ztry)/
            (*this.*probability)(x,y,z));
    if (rand.Rannyu()<=alpha) { //valid even if alpha = 1
        x=xtry;
        y=ytry;
        z=ztry;
        in++;
    }
    tot++;
}

std::vector<double> hydrogen::get_position() const
{
    std::vector<double> pos={x,y,z};
    return pos;
}


double hydrogen::get_radius() const
{
    return std::sqrt(x*x+y*y+z*z);
}


