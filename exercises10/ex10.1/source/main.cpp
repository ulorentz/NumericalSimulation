#include <iostream>
#include <fstream>
#include "travelling_salesman.h"


int main(int argc, char** argv){
    if(argc!=6){
        std::cout << "Program for simulated annealing of a travelling salesman " 
            "problem between 30\ncities.\n\n";
        std::cout << "Usage: " << argv[0] << " <beta> <learning_rate> <n_steps>"
            " <location> <NORM> \n\n";
        std::cout << "Where options [compulsory] are:\n"
                  << "  <beta>\t\tBoltzman weight in the simulated annealing\n"
                  << "  <learning_rate>\tHow much beta is increased in each "
                  "iteration\n"
                  << "  <n_steps>\t\tSteps of algorithm\n"
                  << "  <location>\t\tLocation of the cities, could be "
                  "'circle' or 'squared'\n"
                  << "  <norm>\t\tNorm used to compute the cost, could be 'L1'"
                  " or 'L2'\n\n";
        return 0;
    }

    
    
    AnnealingTSP ann(30,argv[4],argv[5]);
    std::ofstream out(std::string("outputs/")+argv[4]+"_cities_location.dat");
    std::vector<std::array <double , 2 > > locations=ann.getLocations();
    for(auto &it:locations)
        out << it[0] << "\t" << it[1] << std::endl;
    out.close();

    //initial best random path
    SalesmanPath best=ann.getPath();
    std::cout << "Initial best fitness " << argv[4] <<": \t" 
              << best.getFitness() << std::endl;

    out.open(std::string("outputs/")+argv[4]+"_initial_best_path.dat");
    std::vector<unsigned int> best_path = best.getPath();
    for(auto &it : best_path)
        out << it << std::endl;
    out.close();
    double beta=atof(argv[1]);
    double learning=atof(argv[2]);
    unsigned int step = atoi(argv[3]);

    std::vector<double> costs=ann.Anneal(beta, learning, step, true);
    best=ann.getPath();
    std::cout << "Final best fitness "<<argv[4]<<": \t" 
              << best.getFitness() << std::endl;
    out.open(std::string("outputs/")+argv[4]+"_final_best_path.dat");
    best_path = best.getPath();
    for(auto &it : best_path)
        out << it << std::endl;
    out.close();

    out.open(std::string("outputs/")+argv[4]+"_lenght.dat");
    for(auto &it:costs)
        out << it << std::endl;
    out.close();
    std::cout << std::endl;
} 
