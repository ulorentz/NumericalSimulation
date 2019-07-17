#ifndef _MOLECULAR_MC_H_
#define _MOLECULAR_MC_H_
#include <string>
#include <cmath>
#include <fstream>
#include <map>
#include <vector>
#include <iostream>
#include "random.h"             //Random

class MolecularMC{
    public:
        //Default constructor: it needs, in the folder of the executable, the
        //following files:
        //  -"Primes" and "seed.in" in order to initialize the random generator
        //  -"input.dat" that is a file containing the simulation parameter. It
        //    will be read in the following way:
        //
        //    ReadInput >> temp;
        //    ReadInput >> npart;
        //    ReadInput >> rho;
        //    ReadInput >> rcut;
        //    ReadInput >> delta;
        //    ReadInput >> nblk;
        //    ReadInput >> nstep;
        //
        // - "initial_configuration" should be the filename of a file containing 
        //   a configuration valid for the molecular model defined in 
        //   "input.dat". 
        // - istantaneous: if true the program will print istantaneous value of
        //                 p and U, and not just the mean with data blocking.
        MolecularMC(std::string initial_configuration, 
                    bool instantaneous=false);

        //Run the simulation
        void Run(); 
        
        ~MolecularMC();
    private:
        Random rand;
        //vectors of positions
        std::vector<double> x, y, z;


        //maps for energy and pressure
        std::vector<std::string> keys;                  //keys for maps
        std::map<std::string, double> walker;
        std::map<std::string, double> block_average;
        std::map<std::string, double> global_average, global_average2;

        //histogram
        std::vector<double> histo_walker;
        std::vector<double> histo_block_average;
        std::vector<double> histo_global_average, histo_global_average2;

        //blocking parms
        unsigned int blk_norm, accepted, attempted; 

        //Parameters 
        double box, temp,  rho, rcut, delta;
        unsigned int npart, nblk, nstep;
        double beta, vol, vtail, ptail;
        const unsigned int nbins = 100;
        double bin_size;
    
        std::ofstream Gerr, Gave, Epot, Pres;
        std::ofstream ist_pot;
        std::ofstream ist_pres;

        bool istant;
        // METHODS
        void Reset(unsigned int); 
        void Accumulate(); 
        void Averages(unsigned int); 
        void Move(); 
        void ConfFinal(); 
        void ConfXYZ(unsigned int); 
        void Measure(); 
        double Boltzmann(double, double, double, unsigned int); 
        double Pbc(double) const;
        double Error(double,double,unsigned int) const;


};


//Inline methods
//

inline double MolecularMC::Pbc(double r) const{
    return r - box * rint(r/box);
}

inline double MolecularMC::Error(double sum, 
                                 double sum2, 
                                 unsigned int iblk) const{
    return sqrt((sum2/(double)iblk - pow(sum/(double)iblk,2))/(double)iblk);
}
#endif //_MOLECULAR_MC_H_
