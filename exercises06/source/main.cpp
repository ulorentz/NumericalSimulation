#include <iostream>
#include "Ising1D.h"

int main(int argc, char** argv){
    std::string usage =
    "Usage: " + std::string(argv[0]) + " [options] \n"+
    "Where options (compulsory) should be: \n" +
    "\t--one-config     to start from random configuration\n"+
    "\t--restart        to restart simulation from configuration file\n"+
    "\t                 named 'config.final'\n";

    if(argc!=2){
        std::cout << "Program for simulating 1D Ising Model\n\n" <<usage;
        return 0;
    }
    
    if(std::string(argv[1])=="--restart"){
        Ising1D ising("config.final");
        ising.Run();
    }
    else if(std::string(argv[1])=="--one-config"){
        Ising1D ising;
        ising.Run();
    }
    else
        std::cout << "\nUnrecognized option '" << argv[1] << "'\n" <<
            std::endl <<usage;
    return 0;
}
