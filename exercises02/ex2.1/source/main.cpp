#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <numeric>          //std::partial_sum
#include "accept_reject.h"  //AcceptReject
#include "functions.h"      //initialize_seed
#include <cmath>




//Normalization
const double norm=(64*std::sqrt(2)-24*M_PI+std::pow(M_PI,3))/(24.*M_PI);
//Function to integrate
double function(double x){
    return M_PI*std::cos(M_PI*x/2.)/2.;
}

//Distribution
double distribution(double x){
    return std::abs(1.-M_PI*M_PI*x*x/8.);
} 

//g
double g(double x){
    return (1.-M_PI*M_PI*x*x/8.);
}

  

int main(){
    AcceptReject rand;
    //random seed setting
    initialize_seed(rand, "Primes", "seed.in");
    //random number generation
    unsigned int m=10000;       //number of throws
    unsigned int n=100;         //number of blocks
    unsigned int k=m/n;         //throws in each blocks
    std::vector<double> integral(k,0.);

    //PART 1
    for(unsigned int i=0; i<n; ++i){
        for(unsigned int j=0; j<k;++j){
            integral[i]+=function(rand.Rannyu());
        }
        integral[i]/=k;
    }
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

    //PART 2
    
    std::fill(integral.begin(), integral.end(), 0);
    double x;
    for(unsigned int i=0; i<n; ++i){
        for(unsigned int j=0; j<k;++j){
            x=rand.distribution(distribution, 1.);
            //TODO EXPLAIN
            integral[i]+=norm*function(x)/g(x);
        }
        integral[i]/=k;
    }
    
    out.open("integral_importance.txt");
    if(out.fail()){
        std::cerr << "Error opening outfile\n";
        return 1;
    }
    //error
    error=blocking_error(integral);
    for(unsigned int i=0;i<n; ++i)
        out << (i+1)*k <<"\t"<< integral[i] << "\t" << error[i] << std::endl;
    
 
    out.close();
 


 
    return 0;
}
