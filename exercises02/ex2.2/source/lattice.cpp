#include "lattice.h"
#include <cmath>
#include <cstdlib>
#include <iostream>
DiscreteLattice::DiscreteLattice(unsigned int dim) :
    position(3)
{
    for (auto &it : position)
        it=0;
}

DiscreteLattice::~DiscreteLattice(){
}
 
void DiscreteLattice::setPosition(const std::vector<int>& pos){
    if (pos.size()!=position.size())
        throw std::range_error("Given vector has wrong dimension");
    position=pos;
}
 
void DiscreteLattice::move(unsigned int axis, bool direction){
    position[axis-1]+=2*direction-1;
}

int DiscreteLattice::rsquared() const{
    int dist=0;
    for (auto& it: position)
        dist+=it*it;
    return dist;
}

void DiscreteLattice::origin(){
    for(auto &it:position)
        it=0;
}


Dense3DRW::Dense3DRW(double step_size) :
    stepsize(step_size)
{
    for (auto &it : position)
        it=0;
}

Dense3DRW::~Dense3DRW(){
}
 
void Dense3DRW::setPosition(const std::array<double, 3>& pos){

}
 
void Dense3DRW::move(double theta, double phi){
    position[0]+=stepsize*std::sin(theta)*std::cos(phi);
    position[1]+=stepsize*std::sin(theta)*std::sin(phi);
    position[2]+=stepsize*std::cos(theta);
}

double Dense3DRW::rsquared() const{
    double dist=0.;
    for (auto& it: position)
        dist+=it*it;
    return dist;
}

void Dense3DRW::origin(){
    for(auto &it:position)
        it=0.;
}
