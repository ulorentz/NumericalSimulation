#include <iostream>
#include <fstream>
#include "travelling_salesman.h"
#include <fstream>

//#define _CONVERGENCE_
//TODO define in a serious way the convergence!
//
//
void find_convergence(GeneticTSP &);

int main(int argc, char** argv){

/*****************************************************************************/
// Circular part
/*****************************************************************************/
    GeneticTSP circle_tsp(10, 30, "circle", "L2");
    
    //write out the cities location
    std::ofstream out("outputs/circle_cities_location.dat");
    std::vector<std::array <double , 2 > > locations=circle_tsp.getLocations();
    for(auto &it:locations)
        out << it[0] << "\t" << it[1] << std::endl;
    out.close();
    
    //initial best random path
    SalesmanPath best=circle_tsp.getBestPath();
    std::cout << "Initial best fitness circle: \t" 
              << best.getFitness() << std::endl;
    out.open("outputs/circle_initial_best_path.dat");
    std::vector<unsigned int> best_path = best.getPath();
    for(auto &it : best_path)
        out << it << std::endl;
    out.close();

#ifdef _CONVERGENCE_    
    find_convergence(circle_tsp);
#else
    std::ofstream av("av_circle.dat");
    std::ofstream bestout("best_circle.dat");
    for(unsigned int i = 0;i<5000 ; ++i){
        circle_tsp.Evolve();
        bestout<<circle_tsp.getFitness() << std::endl;
        av<<circle_tsp.getAveragedFitness() << std::endl;
    }
    av.close();
    bestout.close();

#endif //_CONVERGENCE_


    best=circle_tsp.getBestPath();
    std::cout << "Final best fitness circle: \t" 
              << best.getFitness() << std::endl;
    out.open("outputs/circle_final_best_path.dat");
    best_path = best.getPath();
    for(auto &it : best_path)
        out << it << std::endl;
    out.close();
    std::cout << std::endl;
/*****************************************************************************/
/* Squared part                                                              */
/*****************************************************************************/
    GeneticTSP squared_tsp(10, 30, "squared", "L2");
    
    //write out the cities location
    out.open("outputs/squared_cities_location.dat");
    locations=squared_tsp.getLocations();
    for(auto &it:locations)
        out << it[0] << "\t" << it[1] << std::endl;
    out.close();
    
    //initial best random path
    best=squared_tsp.getBestPath();
    std::cout << "Initial best squared fitness: \t" 
              << best.getFitness() << std::endl;
    out.open("outputs/squared_initial_best_path.dat");
    best_path = best.getPath();
    for(auto &it : best_path)
        out << it << std::endl;
    out.close();

#ifdef _CONVERGENCE_    
    find_convergence(squared_tsp);
#else
    av.open("av_squared.dat");
    bestout.open("best_squared.dat");
    for(unsigned int i = 0;i<5000 ; ++i){
        squared_tsp.Evolve();
        bestout<<squared_tsp.getFitness() << std::endl;
        av<<squared_tsp.getAveragedFitness() << std::endl;
    }
    av.close();
    bestout.close();


#endif //_CONVERGENCE_


    best=squared_tsp.getBestPath();
    std::cout << "Final best fitness squared: \t" 
              << best.getFitness() << std::endl;
    out.open("outputs/squared_final_best_path.dat");
    best_path = best.getPath();
    for(auto &it : best_path)
        out << it << std::endl;
    out.close();
}



void find_convergence(GeneticTSP &tsp){
    double fit1, fit2;
    fit1=tsp.getBestPath().getFitness();
    for(unsigned int i = 0; ; ++i){
        tsp.Evolve();
        if(i%10000==0){
            fit2=tsp.getBestPath().getFitness();
            
            if(fit1-fit2< 0.00001){
                std::cout << "Convergence reached after " << i << " steps.\n";
                break;
            }
            fit1=fit2;
        }
        if(i==1000000){
            std::cout << "Unable to reach convergence after " << i <<" steps\n";
            exit(0);
        }
    }
}
