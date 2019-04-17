#ifndef _FUNCTIONS_H_
#define _FUNCTIONS_H_
#include "random.h"
#include <iosfwd>
#include <string>
#include <cmath>
#include <vector>


//Function that initializes the seed of the random generator defined in the 
//random.h library. 
//
//Arguments:
//  rand: random generator defined in "random.h" to be initialized
//  primes: filename of primes number
//  seedin: filename of multiple seeds
void 
initialize_seed(Random& rand, std::string primes, std::string seedin);

//Function that computes statistical error with data blocking method. 
//Furthermore it computes incremental statistical average of the input vector
//with blocking method.
//
//Arguments:
//  estimate:   a std::vector where in each position there is the block mean of
//              the quantity to calculate the statistical error of.
//              Note that estimate will be modified in the following way:
//                  estimate[i]+=estimate[i-1];
//                  estimate[i]/=(i+1);
//              -> the function will compute the incremental statistical average
//Returns:
//  vector<double>:     a std::vector where in each position there is the 
//                      statistical uncertanty of the block's estimate. 
//                      For the first block the error is set to zero.
std::vector<double> blocking_error(std::vector<double>& estimate);
#endif //_FUNCTIONS_H_
