#include "accept_reject.h"


AcceptReject::AcceptReject() :
    Random()
{
}
AcceptReject::~AcceptReject(){
}

double AcceptReject::distribution( double (*function)(double), double max){
    double y;
    double x;
    do{
        y=Rannyu()*max;
        x=Rannyu();
    }while(y>function(x));
    return x;
}

