#include "functions.h"
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <numeric>          //std::partial_sum

void 
initialize_seed(Random& rand, std::string primes, std::string seedin){
    int seed[4];
    int p1,p2;
    std::ifstream primes_stream(primes);
    if (primes_stream.fail()){
        std::cerr << "Error opening primes\n";
        exit(1);
    }
    primes_stream >> p1 >> p2;
    primes_stream.close();

    std::ifstream seed_reader(seedin);
    if (seed_reader.fail()){
        std::cerr << "Error opening seed\n";
        exit(1);
    }
    std::string property;
    seed_reader >> property;
    while(!seed_reader.eof()){ //is it usefull to read more than once?
       if( property == "RANDOMSEED" ){
            seed_reader >> seed[0] >> seed[1] >> seed[2] >> seed[3];
            rand.SetRandom(seed,p1,p2);
        }
        seed_reader >> property;
    }
    seed_reader.close();
}

 
std::vector<double> blocking_error(std::vector<double>& estimate)
{
    std::vector<double> error(estimate.size(),0.);
    double last=0., last2=0.;
    last=estimate[0];
    last2=last*last;
    for(unsigned i=1; i<estimate.size(); ++i){
        last+=estimate[i];
        last2+=estimate[i]*estimate[i];
        estimate[i]=last/(i+1);
        error[i]=std::sqrt((last2/(i+1)-estimate[i]*estimate[i])/i);
    }
    return error;
}
