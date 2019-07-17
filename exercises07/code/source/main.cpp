#include <iostream>
#include "MolecularMC.h"

int main(int argc, char** argv){
    std::string usage =
    "Usage: " + std::string(argv[0]) + " [options] \n"+
    "Where options (compulsory) should be: \n" +
    "\t--new-sim        to start from 'config.0' configuration\n"+
    "\t--restart        to restart simulation from 'config.final'\n"+
    "\t                 named 'config.final'\n"+
    "Facoltative options is:\n"+
    "\t--istant         to print istantaneous values\n";

    if(argc==2){
        if(std::string(argv[1])=="--restart"){
            MolecularMC mol("config.final");
            mol.Run();
        }
        else if(std::string(argv[1])=="--new-sim"){
            MolecularMC mol("config.0");
            mol.Run();
        }
        else
            std::cout << "\nUnrecognized option '" << argv[1] << "'\n" <<
                std::endl <<usage;
    }
    else if(argc==3){
        if(std::string(argv[1])=="--istant" or
           std::string(argv[2])=="--istant"){
            if(std::string(argv[1])=="--restart" or
               std::string(argv[2])=="--restart"){
                MolecularMC mol("config.final", true);
                mol.Run();
            }
            else if(std::string(argv[1])=="--new-sim" or
                    std::string(argv[2])=="--new-sim"){
                MolecularMC mol("config.0", true);
                mol.Run();
            }
            else
                std::cout << "\nUnrecognized option '" << argv[1] << "'\n" <<
                    std::endl <<usage;
               
        }
        else
            std::cout << "\nUnrecognized option.\n" <<
                std::endl <<usage;
    }
    else{    
        std::cout << "Program for simulating 1D Ising Model\n\n" <<usage;
        return 0;
    }
    return 0;
}
