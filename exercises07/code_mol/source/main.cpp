#include "molecular_dynamics.h"
#include <iostream>
int main(int argc, char** argv){
    std::string usage =
    "Usage: " + std::string(argv[0]) + " [options] \n"+
    "Where options (compulsory) should be: \n" +
    "\t--one-config     to start from one configuration file named\n"+
    "\t                 'config.0'\n" +
    "\t--restart        to restart simulation from configuration file\n"+
    "\t                 named 'config.final' and 'config.old'\n";

    if(argc!=2){
        std::cout << "Program for computational molecular dynamics\n\n" <<usage;
               return 0;
    }
    
    if(std::string(argv[1])=="--restart"){
        MolecularDynamics molDin("input.dat", "config.final", "config.old");
        molDin.RunSimulation();
    }
    else if(std::string(argv[1])=="--one-config"){
        MolecularDynamics molDin("input.dat", "config.0");
        molDin.RunSimulation();
    }
    else
        std::cout << "\nUnrecognized option '" << argv[1] << "'\n" <<
            std::endl <<usage;
    return 0;
}
