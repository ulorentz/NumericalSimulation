#include <iostream>
#include <fstream>
#include "travelling_salesman.h"
#include <fstream>


//TODO define in a serious way the convergence!
/* run_genetic is a function that simply setup the code for the genetic 
 * algorithm, where:
 *  - topology:  is or "circle" or " squared", and indicates how cities are 
 *    located.
 *  - norm: can be "L2" or "L1"
 *  - n_genes is an integer indicating the number of genes in the population
 *  - n_cities is an integer that specifies the number of cities
 *  The true code is in the library "travelling_salesman", this main simply 
 *  calls functions and classes defined there.
 */




#define _CONVERGENCE_
//void find_convergence(GeneticTSP &);
void find_convergence(GeneticTSP &tsp, std::ofstream& av, std::ofstream& best);

void run_genetic(std::string topology, 
                 std::string norm, 
                 unsigned int n_genes, 
                 unsigned int n_cities)
{
    GeneticTSP tsp(n_genes, n_cities, topology, norm);
    std::ofstream out("outputs/"+topology +"_cities_location.dat");
    std::vector<std::array <double , 2 > > locations;
    locations=tsp.getLocations();
    for(auto &it:locations)
        out << it[0] << "\t" << it[1] << std::endl;
    out.close();
    SalesmanPath best=tsp.getBestPath();
    std::cout << "Initial best fitness" << topology << " " << norm <<": \t"
              << best.getFitness() << std::endl;   
    out.open("outputs/"+topology+"_initial_best_path.dat");
    std::vector<unsigned int> best_path = best.getPath();
    for(auto &it : best_path)
        out << it << std::endl;
    out.close();
    
    
    std::ofstream out_av("outputs/fit_av_"+topology+"_"+norm+".dat");
    std::ofstream out_best("outputs/fit_best_"+topology+"_"+norm+".dat");
    find_convergence(tsp, out_av, out_best);
    out_av.close();
    out_best.close();

    best=tsp.getBestPath(); 
    std::cout << "Final best fitness" <<topology<<" " <<norm <<": \t"
              << best.getFitness() << std::endl;
    out.open("outputs/"+topology+"_"+norm+ "_final_best_path.dat");
    best_path = best.getPath();
    for(auto &it : best_path)
        out << it << std::endl;
    out.close();
    std::cout << std::endl;
}


int main(){
    //Circle 
    //L2
    run_genetic("circle", "L2", 10, 30);
    //L1
    run_genetic("circle", "L1", 10, 30);

    //Square
    //L2
//run_genetic("squared", "L2", 10, 30);
    run_genetic("squared", "L2", 10, 30);
    //L1
    run_genetic("squared", "L1", 10, 30);
}



/*void find_convergence(GeneticTSP &tsp){
#ifdef _CONVERGENCE_    
   double fit1, fit2, fit3,fit4;
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
#else
    for(unsigned int i = 0;i<1000000 ; ++i)
        tsp.Evolve();
#endif //_CONVERGENCE_
}
*/

//evolve the genetic for step time, and save update counter, that counts the 
//total time of evolutions
double evolve(GeneticTSP &tsp, 
              unsigned int step, 
              unsigned int &counter,
              std::ofstream& av, 
              std::ofstream& best){
    double cost1=0., cost2=0.;
    cost1=tsp.getFitness();
    for(unsigned int i=0; i<step; ++i){
        tsp.Evolve();
        av<<tsp.getAveragedFitness() << std::endl;
        best << tsp.getFitness() << std::endl;
    }
    counter+=step;
    cost2=tsp.getFitness();
    return cost1-cost2;
}

void find_convergence(GeneticTSP &tsp, std::ofstream& av, std::ofstream& best){
    //basically I require that convergence is reached if for "tot" times
    //the cost  decreases less than 0.001
    unsigned int counter=0;
    double cost=0.;
    unsigned int sum=0;
    const unsigned int tot=100, breaker=10000000;
    do{
        sum=0;
        for(unsigned int i=0;i<tot; ++i){    
            cost=evolve(tsp, 1000, counter, av, best);
            if(cost<0.001)
                sum+=1;
        }
        if(sum==tot) 
            break;
    }while(counter!=breaker);
    if(counter==breaker)
        std::cout<<"Unable to reach convergence after " << breaker<<" steps\n";
    else
        std::cout <<"Convergence found after " << counter << " steps\n";
}
