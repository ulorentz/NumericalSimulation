#include "random.h"
#include "functions.h"      //initialize_seed
#include <fstream>
#include <iostream>
#include <cmath>
#include <vector>
#include <numeric>           //partial_sum


int main(){
    Random rand;
    initialize_seed(rand, "Primes", "seed.in");
    unsigned int n=(int)1E7;    //total number of point
    unsigned int m=100000;      //point in each block
    unsigned int k=n/m;         //number of blocks
    unsigned int n_hit=0;       //counter of hit
    double D=10., L=4;          //parameters of experiment
    //std::array <double, k> sum_pi, sum_pi2;
    std::vector <double> sum_pi(k);
    double x, len, pi;
    //Random number generation
    std::ofstream out("pi.txt");
    double x1, x2;    
    for(unsigned int j=0; j<k; ++j){
        n_hit=0;
        for(unsigned int i=0; i<m; ++i){
            //computing pi
            x=rand.Rannyu(0, D);    //position of barycenter 
            x1=rand.Rannyu(-1,1);   //x1, x2 are coordinate into rectangle
            x2=rand.Rannyu();
            while(x1*x1+x2*x2>1){   //ensure x1 is distributed as a semicircle
                x1=rand.Rannyu(-1,1);
                x2=rand.Rannyu();
            }
            len=L*x1/std::sqrt(x1*x1+x2*x2);    //normalize x1 to have a cosine
            if ((x-len)<=0. or (x+len)>=D)      //check if lines are hit
                n_hit++;
        }   
        pi=2.*L*m/(n_hit*D);
        sum_pi[j]=pi;
    }
    //computing error and printing
    std::vector<double> error=blocking_error(sum_pi);
    for(unsigned int i=0; i<k; ++i)
        out << (i+1)*m <<"\t"<< sum_pi[i] << "\t" << error[i] << std::endl;
    
    
    
    out.close();
}
