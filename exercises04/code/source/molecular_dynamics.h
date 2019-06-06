#ifndef _MOLECULAR_DYNAMICS_H_
#define _MOLECULAR_DYNAMICS_H_
#include "random.h"
#include <vector>
#include <string>
#include <cmath>
#include <fstream>

class MolecularDynamics{
    public:
        //Input parameters:
        // - 'simParameters' is the filename for input parameters
        // - 'configFile' is the filename for initial molecular configurations
        MolecularDynamics(std::string simParameters, 
                          std::string configFile);
        //Input parameters:
        // - 'simParameters' is the filename for input parameters
        // - 'configFile' is the filename for initial molecular configurations,
        //    usually is the last configuration of a previous simulation.
        // - 'oldConfigFile' is the filename for the molecular configuration 
        //    prior to 'configFile': it is used to extrapolate velocities.
        MolecularDynamics(std::string inputFile, 
                          std::string configFile,
                          std::string oldConfigFile);

        void RunSimulation();
        
        ~MolecularDynamics();

    private:
        /*** Private methods ***/
        void Force(); //compute all forces on all particles
        double Pbc(double) const; 
        void Move();
        void ConfFinal(std::string filename) const;
        void ConfXYZ(unsigned int) const;
        void Measure();
        //compute and print statistical average and errors by blocking method
        void PrintBlocking();
      
        /*** Data members ***/
        //number of particles to simulate
        unsigned int npart;
        //vectors of positions, old positions, velocities and forces. 
        std::vector<double> x, y, z, xold, yold, zold, vx, vy, vz, fx, fy, fz;

        //physical parameters
        double energy,temp,vol,rho,box,rcut;
        //parameters of similation
        double delta;
        unsigned int iprint, nstep;

        Random rand; // random class 

        //parameters, observables
        double  stima_pot, stima_kin, stima_temp, stima_etot, stima_press; 
        //blocking vectors
        std::vector<double> est_pot, est_kin, est_etot, est_temp, est_press;

        //measure time interval and counter
        unsigned int measure_time_interval;
        //size of each block
        unsigned int block_size;
        //actual filling block
        unsigned int imeasure, iblock;
        //output streams
        std::ofstream Epot, Ekin, Etot, Temp, Press;

};

//Algorithm for periodic boundary conditions with side L=box
inline double MolecularDynamics::Pbc(double r) const{  
    return r - box * rint(r/box);
}

#endif //_MOLECULAR_DYNAMICS_H_
