#include "Ising1D.h"
#include <iostream>
#include <cstdlib>
#include <iomanip>

Ising1D::Ising1D(std::string previous_configuration) :
    rnd("Primes", "seed.in"),
    keys{"energy", "capacity", "magnetization", "susceptibility"}
{

    for (auto &it : keys){
        walker[it]=0.;
        block_average[it]=0.;
        global_average[it]=0.; 
        global_average2[it]=0.;
    }
    
    //read "input.dat" and set simulation parameters
    GetInputFromFile();
   
    if(metro){
        Move=&Ising1D::MetropolisMove;
        std::cout << "The program perform Metropolis moves" << std::endl; 
    }
    else{
        Move=&Ising1D::GibbsMove;
        std::cout << "The program perform Gibbs moves" << std::endl;
    }
    
   
    if(previous_configuration==""){ //start from random config
        Ene.open("outputs/ene.dat");   //delete content of file
        Heat.open("outputs/heat.dat");
        Mag.open("outputs/mag.dat"); 
        Chi.open("outputs/chi.dat");  
        //initial configuration
        s.resize(nspin);
        for(auto& it : s){
            if(rnd.Rannyu() >= 0.5) it=1;
            else it = -1;
        }
    }
    else{
        Ene.open("outputs/ene.dat", std::ios::app);   
        Heat.open("outputs/heat.dat", std::ios::app);
        Mag.open("outputs/mag.dat", std::ios::app); 
        Chi.open("outputs/chi.dat", std::ios::app);  
        std::ifstream input_file(previous_configuration);
        if(input_file.fail()){
              std::cout << "Unable to open " << previous_configuration << "\n";
              exit(1);
        }
        short tmp;
        input_file >> tmp;
        while(!input_file.eof()){
            s.push_back(tmp);
            input_file >> tmp;
        }
        input_file.close();
        if (s.size()!=nspin){ 
            std::cout << "Spin configurations found  in '" 
                      << previous_configuration 
                      << "' don't belong to a simulation compatible with "
                      << "parameters defined in 'input.dat'\n";
            exit(0);
        }

    }
    
    //Evaluate energy etc. of the initial configuration
    Measure();

    //Print initial values for the potential energy and virial
    std::cout << "Initial energy = " 
              << walker.at("energy")/(double)nspin << std::endl;

}



void Ising1D::GetInputFromFile()
{
    //Read input informations
    std::ifstream ReadInput("input.dat");
    if(ReadInput.fail()){
        std::cout << "Unable to open 'input.dat'\n";
        exit(1);
    }
 
    std::cout << "Classic 1D Ising model             \n"; 
    std::cout << "Monte Carlo simulation             \n\n"; 
    std::cout << "Nearest neighbour interaction      \n\n"; 
    std::cout << "Boltzmann weight exp(- beta * H ), beta = 1/T \n\n"; 
    std::cout << "The program uses k_B=1 and mu_B=1 units \n";
  
    ReadInput >> temp;
    beta = 1.0/temp;
    std::cout << "Temperature = " << temp << std::endl;
    
    ReadInput >> nspin;
    std::cout << "Number of spins = " << nspin << std::endl;
    
    ReadInput >> J;
    std::cout << "Exchange interaction = " << J << std::endl;
    
    ReadInput >> h;
    std::cout << "External field = " << h << std::endl << std::endl;
    
    //set metropolis or gibbs  
    ReadInput >> metro; // if=1 Metropolis else Gibbs

     ReadInput >> nblk;
    
    ReadInput >> nstep;
    
    std::cout << "Number of blocks = " << nblk << std::endl;
    std::cout << "Number of steps in one block = " << nstep 
              << std::endl << std::endl;
    ReadInput.close();


}

void Ising1D::Measure()
{
    double u = 0.,  m = 0.;
    //cycle over spins
    for (unsigned int i=0; i<nspin; ++i){
        u += -J * s[i] * s[Pbc(i+1)] - 0.5 * h * (s[i] + s[Pbc(i+1)]);
        m += s[i];
    }
    walker.at("energy") = u;
    walker.at("capacity") = u*u;
    walker.at("magnetization") = m;
    walker.at("susceptibility") = m*m;
}

void Ising1D::Reset(unsigned int iblk) //Reset block averages
{
   if(iblk==0)
   {
       for(auto &it : keys)
       {
           global_average.at(it) = 0.;
           global_average2.at(it) = 0.;
       }
   }
   
   for(auto &it : keys)
     block_average.at(it) = 0.;
   
   blk_norm = 0.;
   attempted = 0.;
   accepted = 0.;
}

void Ising1D::Accumulate()
{

    for(auto &it : keys)
        block_average.at(it) = block_average.at(it) + walker.at(it); 
    
    blk_norm++;
}

void Ising1D::BlockAverages(unsigned int iblk) //Print results for current block
{
 
    const int wd=12;
     
    std::cout << "Block number " << iblk << std::endl;
    std::cout << "Acceptance rate " << accepted/attempted << "\n\n";
    
    //Energy
    stima_u = block_average.at("energy")/blk_norm/(double)nspin; 
    global_average.at("energy")  += stima_u;
    global_average2.at("energy") += stima_u*stima_u;
    err_u=Error(global_average.at("energy"),global_average2.at("energy"),iblk);
    
    Ene << std::setw(wd) << iblk <<  std::setw(wd) << stima_u << std::setw(wd) 
        << global_average.at("energy")/(double)iblk << std::setw(wd) 
        << err_u <<std::endl;
    
    //Capacity
    stima_c = block_average.at("capacity")/blk_norm/(double)nspin; 
    stima_c = (stima_c-stima_u*stima_u*(double)nspin)/(temp*temp);
    global_average.at("capacity")  += stima_c;
    global_average2.at("capacity") += stima_c*stima_c;
    err_c=Error(global_average.at("capacity"),
                global_average2.at("capacity"),
                iblk);
    
    Heat<< std::setw(wd) << iblk <<  std::setw(wd) << stima_c << std::setw(wd) 
        << global_average.at("capacity")/(double)iblk << std::setw(wd) 
        << err_c <<std::endl;
    
    //Magnetization
    stima_m = block_average.at("magnetization")/blk_norm/(double)nspin;
    global_average.at("magnetization")  += stima_m;
    global_average2.at("magnetization") += stima_m*stima_m;
    err_m=Error(global_average.at("magnetization"),
                global_average2.at("magnetization"),
                iblk);
    
//    Mag << std::setw(wd) << iblk <<  std::setw(wd) << stima_m << std::setw(wd) 
//        << global_average.at("magnetization")/(double)iblk << std::setw(wd) 
//        <<"\t"<< err_m <<std::endl;
    Mag  << iblk <<  "\t" << stima_m << "\t" 
        << global_average.at("magnetization")/(double)iblk  
        <<"\t"<< err_m <<std::endl;
   
    //Susceptibility
    if(h!=0)
        stima_x=0;
    else
        stima_x =block_average.at("susceptibility")/blk_norm/(double)nspin/temp;
    
    global_average.at("susceptibility")  += stima_x;
    global_average2.at("susceptibility") += stima_x*stima_x;
    err_x=Error(global_average.at("susceptibility"),
                global_average2.at("susceptibility"),
                iblk);
    
    Chi << std::setw(wd) << iblk <<  std::setw(wd) << stima_x << std::setw(wd) 
        << global_average.at("susceptibility")/(double)iblk << std::setw(wd) 
        << err_x <<std::endl;
    std::cout << "----------------------------\n\n"; 
}


void Ising1D::ConfFinal(void)
{
    std::ofstream WriteConf;
    
    std::cout << "Print final configuration to file config.final\n\n";
    WriteConf.open("config.final");
    for (auto &it: s)
        WriteConf << it << std::endl;
    WriteConf.close();
    
    rnd.SaveSeed();
}

void Ising1D::Run()
{
    for(int iblk=1; iblk <= nblk; ++iblk) //Simulation
    {
        Reset(iblk);   //Reset block averages
        for(int istep=1; istep <= nstep; ++istep)
        {
            (*this.*Move)();
            Measure();
            Accumulate(); //Update block averages
        }
        BlockAverages(iblk);   //Print results for current block
    }
    ConfFinal(); //Write final configuration
}

void Ising1D::MetropolisMove()
{
    unsigned int iflip=rnd.Rannyu()*nspin;  
    double energy_diff=2*J*(-1*s[iflip])*(s[Pbc(iflip-1)]+s[Pbc(iflip+1)])-
        2*h*s[iflip]; 
    energy_diff=std::exp(energy_diff/temp);
    attempted++;
    if(energy_diff>=1){
        s[iflip]*=-1;
        accepted++;
    }
    else if(rnd.Rannyu()<energy_diff){
        s[iflip]*=-1;
        accepted++;
    }
}

void Ising1D::GibbsMove()
{
    unsigned int iflip=rnd.Rannyu()*nspin;  
    double e1=Boltzmann(1, iflip); 
    double e2=Boltzmann(-1, iflip); 
    double zeta = std::exp(-beta*e1)+std::exp(-beta*e2);

    double up= std::exp(-beta*e1)/zeta;
    if(up>rnd.Rannyu())
        s[iflip]=1;
    else
        s[iflip]=-1;
}

Ising1D::~Ising1D()
{
    Ene.close();
    Heat.close();
    Mag.close();
    Chi.close();
}
