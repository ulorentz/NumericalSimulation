#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include "functions.h"      //blocking_error
#include <cmath>
#include "hydrogen.h"

void help_message(char* argv[]){
    std::cout <<"Usage: \t" << argv[0] << " [options]\n";
    std::cout << "Options:" << std::endl;
    std::cout <<"\t-h, --help\t\t\tprint out this help message\n" <<std::endl;
    std::cout <<"\t--pos [n-points] [generator]\tprint to file positions and "
       "not mean\n\t\t\t\t\tRadius. If n-points (integer) is given,\n\t\t\t\t\t"
       "only this number of moves is done. \n\t\t\t\t\tDefault is 1E6"
       "\n\t\t\t\t\tGenerator could be 'uniform' or 'gauss',"
       "\n\t\t\t\t\tdefault is uniform.\n";
    exit(0);
}


int main(int argc, char* argv[]){
    //Options setup
    if(argc>4){
        std::cout << "Unrecognized option" << std::endl;
        help_message(argv);
    }
    if (argc == 2 or argc == 3 or argc==4){
        if (std::string(argv[1])=="--help" or std::string(argv[1])=="-h")
            help_message(argv);
        if (std::string(argv[1])!="--pos"){
            std::cout << "Unrecognized option" << std::endl;
            help_message(argv);
        }
    }
    
    unsigned int n=(int)1E6; //number of samples
    std::string generator="uniform";
    //three options
    if (argc == 3){
        char *end;
        n = strtol(argv[2], &end, 10);
        if (*end != '\0') { //unable to convert to int
             if(std::string(argv[2])!="uniform" and 
                std::string(argv[2])!="gauss"){
                std::cout << "Unrecognized generator.\n";
                help_message(argv);
            }
            else
                generator=argv[2];
        }
    }
    if(argc==4){
        char *end;
        n = strtol(argv[2], &end,10);
        if((std::string(argv[3])!="uniform" and 
           std::string(argv[3])!="gauss") or *end !='\0')
        {
            n = strtol(argv[3], &end,10);
            if ((std::string(argv[2])!="uniform" and                                  
                std::string(argv[2])!="gauss") or *end !='\0')
            {
                std::cout << "Unrecognized generator.\n";
                help_message(argv);
            }
            else
                generator=argv[2];
        }
        else
            generator=argv[3];
     }

    //end of options setup
    if(argc==1){
        unsigned int block_size=(int)1E4;
        unsigned int n_blocks=n/block_size;
    
    /* Mean radius computation with uniform random numbers */

        hydrogen hydro(1,0,0,2.4, "100", "uniform");
        std::vector<double> r100(n_blocks,0.), r210(n_blocks,0.);
        for (unsigned int i=0; i< n_blocks; ++i){
            for(unsigned j=0; j<block_size; ++j){
                hydro.sample();
                r100[i]+=hydro.get_radius();
            }
            r100[i]/=block_size;
        }
        std::cout << "Uniform 100 acceptance: " << hydro.acceptance() 
            << std::endl;
    
        hydro.reset_metropolis(1,0,0,5.8,"210", "uniform");
        for (unsigned int i=0; i< n_blocks; ++i){
            for(unsigned j=0; j<block_size; ++j){
                hydro.sample();
                r210[i]+=hydro.get_radius();
            }
            r210[i]/=block_size;
        }
        std::cout << "Uniform 210 acceptance: " << hydro.acceptance() 
            << std::endl;
        //computing errors with blocking method
        std::vector<double> error100=blocking_error(r100);
        std::vector<double> error210=blocking_error(r210);
    
        std::ofstream out("radius100uniform.txt");
        for(unsigned int i=0; i<n_blocks; ++i)
            out << r100[i] << "\t" << error100[i] << std::endl;
        out.close();
    
        out.open("radius210uniform.txt");
        for(unsigned int i=0; i<n_blocks; ++i)
            out << r210[i] << "\t" << error210[i] << std::endl;
        out.close();
        
        /* Mean radius computation with gaussian random numbers */
        std::fill(r100.begin(), r100.end(), 0.); 
        std::fill(r210.begin(), r210.end(), 0.); 
    
        hydro.reset_metropolis(1,0,0,1.5,"100", "gauss");
    
        for (unsigned int i=0; i< n_blocks; ++i){
            for(unsigned j=0; j<block_size; ++j){
                hydro.sample();
                r100[i]+=hydro.get_radius();
            }
            r100[i]/=block_size;
        }
        std::cout << "Gauss 100 acceptance: " << hydro.acceptance() 
            << std::endl;
    
        hydro.reset_metropolis(1,0,0,3.5, "210", "gauss");
        for (unsigned int i=0; i< n_blocks; ++i){
            for(unsigned j=0; j<block_size; ++j){
                hydro.sample();
                r210[i]+=hydro.get_radius();
            }
            r210[i]/=block_size;
        }
        std::cout << "Gauss 210 acceptance: " << hydro.acceptance() 
            << std::endl;
        //computing errors with blocking method
        error100=blocking_error(r100);
        error210=blocking_error(r210);
    
        out.open("radius100gauss.txt");
        for(unsigned int i=0; i<n_blocks; ++i)
            out << r100[i] << "\t" << error100[i] << std::endl;
        out.close();
    
        out.open("radius210gauss.txt");
        for(unsigned int i=0; i<n_blocks; ++i)
            out << r210[i] << "\t" << error210[i] << std::endl;
        out.close();
    }
    else {
        double step=0.;
        if (generator == "uniform")
            step=2.4;
        else
            step=1.5;
        hydrogen hydro(1,0,0,step,"100", generator);

        std::vector<std::vector<double> > pos100(n), pos210(n);
        for(unsigned int i=0; i<n; ++i){
            hydro.sample();
            pos100[i]=hydro.get_position();
        };
        //used for tuning acceptance to 50%
        std::cout << "100 acceptance: " << hydro.acceptance() << std::endl;
    
        if (generator == "uniform")
            step=5.8;
        else
            step=3.5;
        hydro.reset_metropolis(1,0,0,step, "210", generator);

        for(unsigned int i=0; i<n; ++i){ 
            hydro.sample();
            pos210[i]=hydro.get_position();
        }
        std::cout << "210 acceptance: " << hydro.acceptance() << std::endl;
        
        std::ofstream out100("100.txt");
        std::ofstream out210("210.txt");
        for(unsigned int i=0; i<n;++i){
            out100 << pos100[i][0] << "\t" 
                   << pos100[i][1] << "\t" 
                   << pos100[i][2] << "\n";
            out210 << pos210[i][0] << "\t" 
                   << pos210[i][1] << "\t" 
                   << pos210[i][2] << "\n";
        }
        out100.close();
        out210.close();
    }


    return 0;
}
