#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <numeric>          //std::partial_sum
#include <algorithm>        //std::max
#include "functions.h"      //initialize_seed, blocking_error
#include "random.h"
#include <cmath>



int main(){
    Random rand;
    //random seed setting
    initialize_seed(rand, "Primes", "seed.in");
    //random number generation
    unsigned int m=10000;       //number of simulations
    unsigned int n=100;         //number of blocks
    unsigned int d=m/n;         //simulations in each blocks

    //finantial parameter
    double s0=100.;             //asset price at t=0
    double del_t=1.;            //delivery time
    double k=100.;              //strike price
    double r=0.1;               //risk free interest rate
    double sigma=0.25;          //volatility

    //Part 1: directly sampled final price
    std::vector<double> call_price(n,0.), put_price(n,0.); 
    double s; 

    //computing stock prices
    for(unsigned int i=0; i<n; ++i){
        for(unsigned int j=0; j<d; ++j){
            s=s0*std::exp((r-0.5*sigma*sigma)*del_t+sigma*rand.Gauss(0,1));
            call_price[i]+=std::exp(-r*del_t)*std::max(0., s-k);
            put_price[i]+=std::exp(-r*del_t)*std::max(0., k-s);
        }
        call_price[i]/=d;
        put_price[i]/=d;
    }

    //computing errors and results
    std::vector<double> error_call=blocking_error(call_price);
    std::vector<double> error_put=blocking_error(put_price);
    std::ofstream out("option_total.txt");
    if(out.fail()){
        std::cerr << "error opening outfile\n";
        return 1;
    }
    for(unsigned int i=0; i<n;++i)
        out << (i+1)*d <<"\t"<< call_price[i] << "\t" << error_call[i] << "\t" 
            << put_price[i] << "\t" << error_put[i] << std::endl;
    
    out.close(); 

    //Part 2: sampling intermediate prices
    std::fill(call_price.begin(), call_price.end(), 0.); 
    std::fill(put_price.begin(), put_price.end(), 0.); 
    unsigned int steps=100;
    double t_step=del_t/steps;
    for(unsigned int i=0; i<n; ++i){
        for(unsigned int j=0; j<d; ++j){
            s=s0;
            //computng final stock price by steps
            for(unsigned int k=0; k<steps; ++k)
                s=s*std::exp((r-0.5*sigma*sigma)*(t_step)+
                        sigma*rand.Gauss(0,1)*std::sqrt(t_step));

            call_price[i]+=std::exp(-r*del_t)*std::max(0., s-k);
            put_price[i]+=std::exp(-r*del_t)*std::max(0., k-s);
        }
        call_price[i]/=d;
        put_price[i]/=d;
    }
    //computing errors and results
    error_call=blocking_error(call_price);
    error_put=blocking_error(put_price);

    out.open("option_steps.txt");
    if(out.fail()){
        std::cerr << "error opening outfile\n";
        return 1;
    }
    for(unsigned int i=0; i<n; ++i)
       out << (i+1)*d <<"\t"<< call_price[i] << "\t" << error_call[i] << "\t"
           << put_price[i] << "\t" << error_put[i]<< std::endl;
    
    out.close(); 



    return 0;
}
