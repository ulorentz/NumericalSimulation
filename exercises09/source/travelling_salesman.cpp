#include "travelling_salesman.h"
#include <iostream>
#include <cmath>

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
        std::cout <<"'" << location << "' it is not a valid option\n";
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


/************ GeneticTSP ***********/

GeneticTSP::GeneticTSP(unsigned int population_size,
                   unsigned int n_cities, 
                   std::string location,
                   std::string _norm) :
    cities(n_cities, location),
    norm(_norm),
    rand("Primes", "seed.in")
{
    for(unsigned int i=0; i<population_size; ++i)
        paths.push_back(SalesmanPath(n_cities));
    for(auto & it : paths)
        it.ComputeFitness(cities, norm);
}

unsigned int GeneticTSP::selection(){
    return paths.size()*std::pow(rand.Rannyu(),2);
}

void GeneticTSP::Evolve(){
    if(rand.Rannyu()<0.99) { //crossover
        //the two parents
        unsigned int sel1=selection();
        unsigned int sel2=selection();
        std::vector<unsigned int> p1=paths[sel1].getPath();
        std::vector<unsigned int> p2=paths[sel2].getPath();
        //decide at wich point break the gene
        unsigned int break_idx=int(rand.Rannyu()*p1.size());
        std::vector<unsigned int> index1;
        std::vector<unsigned int> index2;
        //vector one, let's find in the second vector the missing part of the 
        //first one and do the same for the second one
        std::vector<unsigned int>::iterator it;
        for(unsigned int i=break_idx; i<p1.size(); ++i){
            it=std::find(p2.begin(), p2.end(), p1[i]);
            index2.push_back(it-p2.begin());
            it=std::find(p1.begin(), p1.end(), p2[i]);
            index1.push_back(it-p1.begin());
        }
        //sort the indeces
        std::sort(index1.begin(), index1.end());
        std::sort(index2.begin(), index2.end());

        //fill the ending of the vector with the numbers in the order found.
        std::vector<unsigned int> new_p1, new_p2;
        new_p1=p1;
        new_p2=p2;
        for(unsigned int i=break_idx; i<p1.size(); ++i){
            new_p1[i]=p2[index2[i-break_idx]]; 
            new_p2[i]=p1[index1[i-break_idx]]; 
        }

        //now create the childrens with the given paths
        SalesmanPath path1(new_p1);
        SalesmanPath path2(new_p2);
        path1.ComputeFitness(cities, norm);
        path2.ComputeFitness(cities, norm);
        //add in the end of the queue the two childrens
        paths[paths.size()-2]=path1;
        paths[paths.size()-1]=path2;
    }
    else{ //mutation
        unsigned int sel=selection();
        SalesmanPath p=paths[sel];
        p.Mutation();
        p.ComputeFitness(cities, norm);
        paths[paths.size()-1]=p;
    }
    std::sort(paths.begin(), paths.end());
}


double GeneticTSP::getAveragedFitness() const{
    double av=0.;
    for(unsigned int i=0; i<paths.size()/2; ++i)
        av+=paths[i].getFitness();
    return av/int(paths.size()/2);
}

