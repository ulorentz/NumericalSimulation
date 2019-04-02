#include <fstream>
#include <iostream>
#include "functions.h"

int main(){
    std::ofstream out("rr.txt");
    Random rand;
    
    initialize_seed(rand, "Primes", "seed.in");
    for(unsigned int i=0; i<10000; ++i){
        out <<taylor_p(rand) << std::endl;
    }
    out.close();
}

