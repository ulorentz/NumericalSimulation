#include "distributions.h"
#include "functions.h"
#include <vector>
#include <fstream>
#include <array>
#include <iostream>
#include <numeric>          //std::accumulate
int main(){
    RandomDistributions rand;
    initialize_seed(rand, "Primes", "seed.in");
    std::array<const unsigned int,3> m={2,10,100};
    unsigned int n=10000;
    std::vector <double> gauss(n*m.back()), exp(n*m.back()), cauchy(n*m.back());
    //Random number generation
    for(unsigned int i=0; i< n*m.back(); ++i){
         gauss[i]=rand.Gauss(0.,1.);
         exp[i]=rand.Exponential(1.);
         cauchy[i]=rand.Cauchy(0.,1.);
    }

    std::ofstream out("distributions.txt");
    if(out.fail()){
       std::cerr << "Error opening out file\n";
       return 1;
    }

    for (unsigned int i=0; i<n; ++i){
        //print sampling from distributions
        out << gauss[i] << "\t" << exp[i] << "\t" << cauchy[i];
        //sum and print sampling with M=2,10,100
        for(auto& it:m)
            out << "\t" <<
            std::accumulate(gauss.begin()+i*it, gauss.begin()+(i+1)*it,0.)/it
            << "\t" <<
            std::accumulate(exp.begin()+i*it, exp.begin()+(i+1)*it,0.)/it
            << "\t" <<
            std::accumulate(cauchy.begin()+i*it, cauchy.begin()+(i+1)*it,0.)/it;
        out << std::endl;
    }
    out.close();

}
