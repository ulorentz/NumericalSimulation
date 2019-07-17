#ifndef _TRAVELLING_SALESMAN_H_
#define _TRAVELLING_SALESMAN_H_
#include <vector>
#include <string>
#include <array>
#include <algorithm>
#include "random.h"
//*****************************************************************************/
// In this header three classes are defined:
// - AnnealingTSP,
// - SalesmanPath, 
// - CitiesLocation.
// 
// AnnealingTSP is the class that performs the actual annealing algorithm. To do 
// so, it uses two data types: CitiesLocation, that represents the topology of 
// the cities and SalesmanPath, that represents the path performed by the 
// travelling salesman. 
//*****************************************************************************/

class CitiesLocation{
    public:
        // - n_cities: specifies how many cities there should be
        // - location: options are "circle" and "squared". The first case all 
        //             the cities are located randomly on a circle. In the 
        //             second case all cities are located randomly inside a 
        //             square.
        CitiesLocation(unsigned int n_cities, std::string location="circle");

        //Returns L1 distance between cities a and b, with a starting from 1 
        //and b ending to "n_cities"
        double getL1(unsigned int a, unsigned int b) const;

        //Returns L2 distance between cities a and b, with a starting from 1 
        //and b ending to "n_cities"
        double getL2(unsigned int a, unsigned int b) const;
        
        //Returns positions. Is a vector of couples of coordinate.
        std::vector<std::array <double , 2 > > getLocations() const;
    private:
        //a vector containing ordered cities location
        std::vector<std::array < double, 2> > positions;

        //matrix of distances, computed with L2 and with L1 norm
        std::vector< std::vector <double> > l1, l2;
        //a random object used to decide positions of cities
        Random rand;

};


class SalesmanPath{
    //Note: due to the default indexing I prefer to have the path starting from
    //0, so the completely ordered path would be [0,1,2,...,n_cities-1]
    public:
        //Fill the path with a random valid path, i.e.: [0,5,3,...,n_cities]
        //To be valid all cities must be visited and no city can be explored 
        //more than once.
        //All this check are performed in the private member "check_valid()".
        SalesmanPath(unsigned int n_cities);

        //copy a given path into the salesman path
        SalesmanPath(const std::vector <unsigned int>&);
        ~SalesmanPath();

        //Erase actual path and fill with a new random one
        void RandomPath();

        //compute fitness of path
        // - loc is an object of the class CitiesLocation that specifies where
        //   cities are located
        // - norm is the type of norm to be applied. Options are "L2" and "L1"
        //fitness is stored in a local variable and can be accessed with method
        //"GetFitness"
        void ComputeFitness(const CitiesLocation & loc, std::string norm="L2"); 
        
        //used to access to fitness
        double getFitness() const;

        //used to access the path
        std::vector<unsigned int> getPath() const;


        //perform a mutation on a random allele
        //the type of mutation is chosen at random
        void Mutation();

        //overload of relation operators, relations are based on fitness value
        friend bool operator<(const SalesmanPath &, const SalesmanPath &);
        friend bool operator==(const SalesmanPath &, const SalesmanPath &);
        
        //returns true if the path fulfill the constraints of the salesman 
        //problem, false elsewhere. 
        bool CheckValid() const;
    private:
        std::vector<unsigned int> path;
        double fitness;
        //static random generator  
        static Random rand;

        //methods       
        //
        void single_permutation();

        //all cities are shifted of m position, with m chosen at random
        void shift();

        //m cities (m<=n_cities/2) are shifted at the end of vector
        void partial_shift();

        //a permutation of m contiguous cities, m<n_cities/2 is chosen at random
        void multi_permutation(); 
        
        //inversion of j contiguous cities starting at position k, both numbers
        //are chosen at random
        void inversion();

};


class AnnealingTSP{
    public:
        //Algorithm for simulated annealing of the travelling salesman problem 
        //Arguments:
        // - n_cities: the number of cities of the path
        // - location: how the cities are located, options are "circle" and
        //             "squared"
        // - norm: norm to be applied for computing fitness. 
        //         Options are "L1" and "L2".
        AnnealingTSP(unsigned int n_cities, 
                     std::string location="circle",
                     std::string norm="L2");

        //Get locations of the cities
        std::vector<std::array <double , 2 > > getLocations() const;


        //returns the best path
        SalesmanPath getPath();
        
        //returns the best fitness
        double getFitness() const;


        //simulate the annealing
        // - beta: starting beta in the "hamiltonian" Exp(-b Cost)
        // - learning rate: how much beta is increased on each step
        // - steps: number of steps
        // - return_costs: if true returns a vector containing cost of the path
        //                 for each iteration. 
        //returns:
        // std::vector<double> containing cost of each iteration, if 
        //                     "return_costs" has been set to true.
        std::vector<double> Anneal(double beta, 
                                   double learning_rate, 
                                   unsigned int steps,
                                   bool return_costs=false);
        
        //reset the path with a random path
        void Reset();
    private:
        CitiesLocation cities;
        std::string norm;
        SalesmanPath path;
        Random rand;
};





////////////////////
// Inline methods //
////////////////////

/******* CitiesLocation ******/
inline double CitiesLocation::getL1(unsigned int a, unsigned int b) const{
    return l1[a][b];
}

inline double CitiesLocation::getL2(unsigned int a, unsigned int b) const{
    return l2[a][b];
}

inline std::vector<std::array <double , 2 > > 
CitiesLocation::getLocations() const{
    return positions;

}

/***** SalesmanPath *****/

inline double SalesmanPath::getFitness() const{
    return fitness;
}

inline std::vector<unsigned int> SalesmanPath::getPath() const{
    return path;
}

inline bool operator<(const SalesmanPath & a, const SalesmanPath & b){
    return a.getFitness()<b.getFitness();
}

inline bool operator==(const SalesmanPath &a, const SalesmanPath &b){
    return a.getFitness()==b.getFitness();
}


/****** AnnealingTSP ******/
inline SalesmanPath AnnealingTSP::getPath(){
        return path;
}

inline std::vector<std::array <double , 2 > > 
AnnealingTSP::getLocations() const{
    return cities.getLocations();
}

inline double AnnealingTSP::getFitness() const{
    return path.getFitness();
}

#endif // _TRAVELLING_SALESMAN_H_
