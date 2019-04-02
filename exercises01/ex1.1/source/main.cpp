#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <numeric>      //std::accumulate
#include "random.h"
#include "functions.h"
#include <cmath>

int main(){
    Random rand;
    //random seed setting
    initialize_seed(rand, "Primes", "seed.in");

    //random number generation
    unsigned int m=100000;      //number of throws
    unsigned int n=100;         //number of blocks
    unsigned int l=m/n;         //throws in each blocks

    /*
     * Note about the choose of coding style:
     * With std::vector and STL functions, instead of dynamic array and multiple
     * for statements, when compiled with -O3, the program is approximatly twice 
     * faster than the pure dynamic array version (moreover the code is quite 
     * compact). 
     */
    

    std::vector<double> r(m);   //random numbers
    for (unsigned int i=0; i<m; ++i)
        r[i]=rand.Rannyu();
 
    //PART 1
    std::vector<double> avg(n), avg2(n);
    for (unsigned int i=0; i<n; ++i){
        //sum of random number for each block
        avg[i]=std::accumulate(r.begin()+i*l, r.begin()+(i+1)*l, 0.)/l;
        avg2[i]=avg[i]*avg[i];
    }

    std::vector<double> progress_sum(n), progress_sum2(n), progress_err(n);
    progress_err[0]=0.; //zero position of stat error set to zero
    for (unsigned int i=0; i<n; ++i){
        //progressive sum
        progress_sum[i]=std::accumulate(avg.begin(), avg.begin()+(i+1),
                0.)/(i+1.);
        progress_sum2[i]=std::accumulate(avg2.begin(), avg2.begin()+(i+1),
                0.)/(i+1.);
        if(i!=0) //statistical error
            progress_err[i]=sqrt((progress_sum2[i]-progress_sum[i]*
                        progress_sum[i])/(double)i);
    }
    //printing results
    std::ofstream out("sampling1.txt");
    if(out.fail()){
        std::cerr << "Error opening output file\n";
        return 2;
    }
    for(unsigned int i=0; i<n; ++i)
        out << i*l << "\t" << progress_sum[i] << "\t" << progress_err[i] <<"\n";
    out.close();
 
    //PART 2
    //sum and sum2
    for (unsigned int i=0; i<n; ++i){
        //sum of random number for each block with formula (r-1/2)**2
        //the formula is in lambda expr
        avg[i]=std::accumulate(r.begin()+i*l, r.begin()+(i+1)*l, 0., 
                [](double x, double y){return x+(y-0.5)*(y-0.5);})/l;
        avg2[i]=avg[i]*avg[i];
    }

    progress_err[0]=0.; //zero position of stat error set to zero
    for (unsigned int i=0; i<n; ++i){
        //progressive sum
        progress_sum[i]=std::accumulate(avg.begin(), avg.begin()+(i+1),
                0.)/(i+1.);
        progress_sum2[i]=std::accumulate(avg2.begin(), avg2.begin()+(i+1),
                0.)/(i+1.);
        if(i!=0) //statistical error
            progress_err[i]=sqrt((progress_sum2[i]-progress_sum[i]*
                        progress_sum[i])/(double)i);
    }

    //printing results
    out.open("sampling2.txt");
    if(out.fail()){
        std::cerr << "Error opening output file\n";
        return 2;
    }
    for(unsigned int i=0; i<n; ++i)
        out << i*l << "\t" << progress_sum[i] << "\t" << progress_err[i] <<"\n";
    out.close();
  
    //PART 3
    n=10000;
    m=100;
    r.resize(m*n);
    for (unsigned int i=0; i<m*n; ++i)
        r[i]=rand.Rannyu();
    unsigned int k=100;

    std::vector<unsigned int> interval(m);
    std::vector <double> chi(k);
    for (unsigned int i=0; i<k; ++i){
        std::fill(interval.begin(), interval.end(), 0);
        for(unsigned int j=n*i; j<n*(i+1); ++j)
            interval[(unsigned int)(r[j]*m)]++;
        chi[i]=std::accumulate(interval.begin(), interval.end(), 0.,
                [=](double x, double y){return x+(y-n/m)*(y-n/m)*m/n;});
    }
    //printing results
    out.open("chi.txt");
    if(out.fail()){
        std::cerr << "Error opening output file\n";
        return 2;
    }
    for(auto& it : chi)
        out << it << std::endl;
    out.close();

    return 0;
}
