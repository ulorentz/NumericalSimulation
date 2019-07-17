#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include "functions.h"      //blocking_error
#include "variationalMC.h"

int main(){
    Particle1D Part(0.8888,0.6648,0.,2.7);
   // return 1;
    unsigned int m=10000;       //number of throws
    unsigned int n=50;         //number of blocks
    unsigned int k=m/n;         //throws in each blocks
    std::vector<double> integral(k,0.);
    std::vector<double> psi2;
    //PART 1
    double sample;
    for(unsigned int i=0; i<n; ++i){
        for(unsigned int j=0; j<k;++j){
            sample=Part.Sample();
            integral[i]+=Part.Energy(sample);
            psi2.push_back(sample);
        }
        integral[i]/=k;
    }
    std::cout << "Metropolis sampling\n";
    std::cout << "Acceptance rate: " << Part.Acceptance() << std::endl;
    std::ofstream out("integral.txt");
    if(out.fail()){
        std::cerr << "Error opening outfile\n";
        return 1;
    }
    //compute error and print
    std::vector<double> error=blocking_error(integral);
    for(unsigned int i=0;i<n; ++i)
        out << (i+1)*k <<"\t"<< integral[i] << "\t" << error[i] << std::endl;
    out.close();
    out.open("psi2.txt");
    if(out.fail()){
        std::cerr << "Error opening outfile\n";
        return 1;
    }
    for(auto &it : psi2)
        out << it << std::endl;
    out.close();

    double sigma=0, 
           mu=0,
           var_sigma=1.,
           var_mu=1.,
           start_x=0.,
           jump=3.;
    unsigned int n_samples=10,
                 n_integral=10000,
                 iterations=3;
    Part.FindOptimum(mu,
                     var_mu, 
                     sigma,
                     var_sigma, 
                     start_x, 
                     jump,
                     n_samples,
                     n_integral,
                     iterations);
    Part.GetOptimalParameters(mu, sigma);
    std::cout << "Mu: " << mu << "\nSigma: " << sigma << std::endl;
    std::cout << "Acceptance rate: " << Part.Acceptance() << std::endl;

}

