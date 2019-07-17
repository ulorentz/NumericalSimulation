#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include "mpi.h"
#include "travelling_salesman.h"

#define N_CITIES 30
int main(int argc, char **argv){
     MPI::Init(argc, argv);
     int size = MPI::COMM_WORLD.Get_size();
     int rank = MPI::COMM_WORLD.Get_rank();

     if(argc!=6){
        if(rank==0){
        std::cout << "Argc " << argc << std::endl;
        std::cout << "Program for simulated annealing of a travelling salesman " 
            "problem between 30\ncities.\n\n";
        std::cout << "Usage: " << argv[0] << " <beta> <learning_rate> <n_steps>"
            " <location> <norm> \n\n";
        std::cout << "Where options [compulsory] are:\n"
                  << "  <beta>\t\tBoltzman weight in the simulated annealing\n"
                  << "  <learning_rate>\tHow much beta is increased in each "
                  "iteration\n"
                  << "  <n_steps>\t\tSteps of algorithm\n"
                  << "  <location>\t\tLocation of the cities, could be "
                  "'circle' or 'squared'\n"
                  << "  <norm>\t\tNorm used to compute the cost, could be 'L1'"
                  " or 'L2'\n\n";
        }

        MPI::Finalize();
        return 0;
    }


    //to store the cost of each thread
    double *costs=0;
    if (rank==0)
       costs = new double[size];

    //a unique location of the cities that will be used by all the threads
    double x[N_CITIES];
    double y[N_CITIES];

    //create the location of the cities that will be shared between all threads
    if(rank==0){
        CitiesLocation cities(N_CITIES, argv[4]);
        std::vector<double> tmp=cities.getX();
        for(unsigned int i=0;i<N_CITIES;++i)
            x[i]=tmp[i];
        
        tmp=cities.getY();
        for(unsigned int i=0;i<N_CITIES;++i)
            y[i]=tmp[i];
        //writing cities location
        std::ofstream out(std::string("outputs/")+argv[4]+
                          "_cities_location.dat");
        for(unsigned int i=0; i<N_CITIES; ++i)
            out << x[i] << "\t" << y[i] << std::endl;
        out.close();
    }   
    //before proceding thread 0 should finishes reading    
    MPI_Barrier(MPI::COMM_WORLD);
    //send to all the x and y
    MPI_Bcast(x,N_CITIES,MPI_REAL8, 0, MPI::COMM_WORLD);
    MPI_Bcast(y,N_CITIES,MPI_REAL8, 0, MPI::COMM_WORLD);

    CitiesLocation cities(x,y, N_CITIES);
    AnnealingTSP ann(cities,argv[5], rank);
    
    // add some noise to the random search?
    Random rand("Primes", "seed.in", rank);
    double beta=atof(argv[1]);
    beta=rand.Gauss(beta, (beta/10)*(rank+1));
    double learning=atof(argv[2]);
    learning=rand.Gauss(learning, (learning/10)*(rank+1));
    unsigned int step = atoi(argv[3]);

    ann.Anneal(beta, learning, step, false);

    //communicate costs to the thread 0
    double single_cost=ann.getFitness();
    std::cout << "My id: " << rank << "\t cost: " << single_cost<<std::endl;
    MPI_Gather(&single_cost,1,MPI_REAL8,costs,1,MPI_REAL8,0,MPI::COMM_WORLD);

    int best_rank;
    //find the best TSP path
    if(rank==0){
        best_rank=std::min_element(costs, costs+size)-costs;
    }
    MPI_Bcast(&best_rank, 1, MPI_INTEGER, 0, MPI::COMM_WORLD); 

    //write to disk
    if(rank==best_rank){
        std::cout << "Best thread: " << rank << std::endl;
        SalesmanPath best=ann.getPath();
        std::cout << "Final best fitness "<<argv[4]<<": \t" 
                  << best.getFitness() << std::endl;
        std::ofstream out(std::string("outputs/")+argv[4]+
                          "_final_best_path.dat");
        std::vector<unsigned int> best_path = best.getPath();
        for(auto &it : best_path)
            out << it << std::endl;
        out.close();
    }
   
    if(rank==0)
        delete[] costs;
    MPI::Finalize();


}

