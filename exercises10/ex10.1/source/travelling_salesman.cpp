#include "travelling_salesman.h"
#include <iostream>
#include <cmath>
#include <cassert>


/*********** CitiesLocation ***********/
CitiesLocation::CitiesLocation(unsigned int n_cities, 
                               std::string location) :
    positions(n_cities),
    l1(n_cities, std::vector<double>(n_cities)),
    l2(n_cities, std::vector<double>(n_cities)),
    rand("Primes", "seed.in")
{
    
    if(location == "circle"){//i am supposing the circ has radius = 1
        double theta;
        for (unsigned int i =0; i<n_cities; ++i){
            theta=rand.Rannyu()*M_PI*2.;
            positions[i][0]=std::cos(theta);
            positions[i][1]=std::sin(theta);
        }
    }
    else if(location == "squared"){
        for (unsigned int i =0; i<n_cities; ++i){
            positions[i][0]=rand.Rannyu();
            positions[i][1]=rand.Rannyu();
        }
    }
    else{
        std::cout <<"'" << location << "' it is not a valid location option\n";
        exit(0);
    }

    //l1 distances and l2 distances
    double dist=0.;
    std::array<double, 2> diff;
    for(unsigned int i=0; i< n_cities-1; ++i){ //due to the trick used for the 
                                               //simmetry of the matrix the last
                                               //element is not needed
        for(unsigned int j=i+1; j<n_cities; ++j){
            diff[0]=positions[i][0]-positions[j][0];
            diff[1]=positions[i][1]-positions[j][1];
            dist=diff[0]*diff[0]+diff[1]*diff[1];
            //since the matrix of distances is simmetric
            l2[i][j]=dist;
            l2[j][i]=dist;
            dist=std::sqrt(dist);
            l1[i][j]=dist;
            l1[j][i]=dist;
        }
    }
    //is it useful?
    for(unsigned int i=0; i<n_cities; ++i){
        l1[i][i]=0.;
        l2[i][i]=0.;
    }

}

/************ SalesmanPath *************/
//static member
Random SalesmanPath::rand("Primes", "seed.in");

SalesmanPath::SalesmanPath(unsigned int n_cities) :
    path(n_cities)
{   //an ordered path
    for(unsigned int i=0; i<n_cities; ++i)
        path[i]=i;
    std::random_shuffle(path.begin(), path.end()); //a generic random gene
    fitness=0.; 
}

SalesmanPath::SalesmanPath(const std::vector <unsigned int>& _path) :
    path(_path)
{
    fitness=0.; 
}

SalesmanPath::~SalesmanPath()
{
}

void SalesmanPath::RandomPath(){
    for(unsigned int i=0; i<path.size(); ++i)
        path[i]=i;
    std::random_shuffle(path.begin(), path.end()); //a generic random gene
    fitness=0.; 
}

void SalesmanPath::ComputeFitness(const CitiesLocation & loc, std::string norm){
    fitness=0.;
    if (norm =="L2"){
        for(unsigned int i=0; i<path.size()-1; ++i)
            fitness+=loc.getL2(path[i], path[i+1]);
        fitness+=loc.getL2(path[path.size()-1], path[0]); //return to first city
    }
    else if (norm == "L1") {
        for(unsigned int i=0; i<path.size()-1; ++i)
            fitness+=loc.getL1(path[i], path[i+1]);
        fitness+=loc.getL1(path[path.size()-1], path[0]); //return to first city
    }
}

void SalesmanPath::Mutation(){
    unsigned int index=int(rand.Rannyu()*5);
    switch(index){
        case 0:
            single_permutation();
            break;
        case 1:
            shift();
            break;
        case 2:
            partial_shift();
            break;
        case 3:
            multi_permutation();
            break;
        case 4:
            inversion();
            break;
    };
}

bool SalesmanPath::CheckValid() const{
    unsigned int total=0, sum=0; 
    for(unsigned int i=0; i<path.size()-1; ++i){
        total+=i;
        sum+=path[i];
        if(std::find(path.begin()+i+1, path.end(), path[i])!=path.end())
            return false; //a duplicate has been found
    }
    total+=path.size()-1;
    sum+=*(path.end()-1);
    if(sum!=total) //not all cities have been visited
        return false;

    return true;
}

void SalesmanPath::single_permutation(){
    unsigned int i, j, tmp;
    i=int(rand.Rannyu()*path.size());
    //TODO is it necessary, or correct at least?
    do{
        j=int(rand.Rannyu()*path.size()); //choose a different index
    }while(i==j);
    tmp=path[i];
    path[i]=path[j];
    path[j]=tmp;
}

void SalesmanPath::shift(){
    //how many position should the vector be shifted, from 1 to path.size-1
    unsigned int pos=int(rand.Rannyu()*(path.size()-1))+1; 
    std::rotate ( path.begin(), path.begin()+pos, path.end() ); 
}

void SalesmanPath::partial_shift(){
    //shift in the end of vector
    unsigned int m=rand.Rannyu()*(path.size()/2); //it's a my own decision to 
                                                  //shift at most size/2 elems
    unsigned int first=rand.Rannyu()*(path.size()/2); //from where
    
    std::vector<unsigned int> tmp(m);
    for(unsigned i=0; i<m; ++i)
        tmp[i]=path[i+first];
    for(unsigned i=first+m; i<path.size(); ++i)
        path[i-m]=path[i];
    for(unsigned i=0; i<m; ++i)
        path[path.size()-m-1+i]=tmp[i];
}


void SalesmanPath::multi_permutation(){
    unsigned int m=rand.Rannyu()*(path.size()/2);
    unsigned int pos1=rand.Rannyu()*(path.size()/2-m);
    unsigned int pos2=rand.Rannyu()*(path.size()/2-m)+path.size()/2;
    unsigned int tmp;
    for(unsigned int i=0; i<m; ++i){
        tmp=path[i+pos1];
        path[i+pos1]=path[i+pos2];
        path[i+pos2]=tmp;
    }
    
}
        
void SalesmanPath::inversion(){
    //how many to reverse
    unsigned int n=int(rand.Rannyu()*(path.size()-1))+2; //from 2 to n_cities  
    unsigned int pos=int(rand.Rannyu()*(path.size()-n));
    std::reverse(path.begin()+pos, path.begin()+pos+n);
}

/********** AnnealingTSP ***********/
AnnealingTSP::AnnealingTSP(unsigned int n_cities, 
                           std::string location,
                           std::string _norm) :
    cities(n_cities, location),
    norm(_norm),
    path(n_cities),
    rand("Primes", "seed.in")
{
    if(norm != "L1" and norm != "L2"){
        std::cout << "'"<<norm <<"' is not a valid norm.\n";
        exit(0);
    }
    path.ComputeFitness(cities, norm);
}


void AnnealingTSP::Reset(){
    path.RandomPath();
    path.ComputeFitness(cities, norm);
}

std::vector<double> AnnealingTSP::Anneal(double beta, 
                                         double learning_rate, 
                                         unsigned int steps,
                                         bool return_costs)
{
    std::vector <double> costs;
    if(return_costs)
        costs.resize(steps);

    double probability=0.;
    for(unsigned int i=0; i<steps; ++i){
        SalesmanPath old_path=path;
        path.Mutation();    //perform a mutation to accept or not in 
                            //the metropolis
        path.ComputeFitness(cities, norm);
        probability=std::exp(-beta*(path.getFitness()-old_path.getFitness()));
        if(rand.Rannyu()>probability) //not accept!
            path=old_path;
        if(return_costs)
            costs[i]=path.getFitness();

        beta+=learning_rate;
    }
    return costs;
}


