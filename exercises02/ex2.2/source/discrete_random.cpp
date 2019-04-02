#include "discrete_random.h"

DiscreteRandom::DiscreteRandom() :
    Random()
{
}

DiscreteRandom::~DiscreteRandom(){
}

int DiscreteRandom::uniform(int a, int b){
   return (int)Rannyu(a, b+1);
}
