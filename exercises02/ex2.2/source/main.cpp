#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <numeric>          //std::partial_sum
#include "functions.h"      //initialize_seed
#include "discrete_random.h"
#include "lattice.h"
#include <cmath>
  

int main(){
    DiscreteRandom rand;
    //random seed setting
    initialize_seed(rand, "Primes", "seed.in");
    //random number generation
    unsigned int m=10000;       //number of simulations
    unsigned int n=100;         //number of blocks
    unsigned int d=m/n;         //simulations in each blocks
    unsigned int moves=100;
    //PART 1
    DiscreteLattice lattice;
    std::vector< std::vector<double> > dist(moves, std::vector<double> (n, 0.));
    std::vector< double> error;
    for(unsigned int j=0; j<n; ++j){
        for(unsigned int k =0; k<d; ++k){
            for(unsigned int i=0; i<moves; ++i){
                lattice.move(rand.uniform(1,3), rand.uniform(0,1));
                dist[i][j]+=lattice.rsquared();
            }
        lattice.origin();
        }
        for(unsigned int i=0; i< moves; ++i)
            dist[i][j]/=d;
    }

    //outfile
    std::ofstream out("discrete_lattice.txt");
    if(out.fail()){
        std::cerr << "Error opening outfile\n";
        return 1;
    }

    double sqrtmodr=0.;
    for(unsigned int i=0; i<moves; ++i){
        error=blocking_error(dist[i]);
        //error propagation
        sqrtmodr=std::sqrt(dist[i][n-1]);
        out << (i+1) <<"\t"<< sqrtmodr << "\t" << 0.5*error[error.size()-1]/
            sqrtmodr << "\n";
    }
 
    out.close();
 
    //PART2
    Dense3DRW random_walk;
    //set to zero both matrices
    for(auto &it : dist)
        std::fill(it.begin(), it.end(), 0.);
    
    for(unsigned int j=0; j<n; ++j){
        for(unsigned int k =0; k<d; ++k){
            for(unsigned int i=0; i<moves; ++i){
                random_walk.move(rand.Rannyu(0,M_PI), rand.Rannyu(0,2*M_PI));
                dist[i][j]+=random_walk.rsquared();
            }
        random_walk.origin();
        }
       for(unsigned int i=0; i< moves; ++i)
            dist[i][j]/=d;
    }

    out.open("dense_lattice.txt");
    if(out.fail()){
        std::cerr << "Error opening outfile\n";
        return 1;
    }

    for(unsigned int i=0; i<moves; ++i){
        error=blocking_error(dist[i]);
        //error propagation
        sqrtmodr=std::sqrt(dist[i][n-1]);
        out << (i+1) <<"\t"<< sqrtmodr << "\t" << 0.5*error[error.size()-1]/
            sqrtmodr << "\n";
    }
 
    out.close();
    
    return 0;
}
