#include "MolecularMC.h"
#include <iomanip>

MolecularMC::MolecularMC(std::string initial_configuration,
                         bool istantaneous) : 
    rand("Primes", "seed.in"),
    keys{"energy", "virial"},
    istant(istantaneous)
{


    std::ifstream ReadInput("input.dat");
    if(ReadInput.fail()){
        std::cerr << "No valid 'input.dat' parameters file found. Aborting.\n";
        exit(1);
    }
    std::ifstream ReadConf(initial_configuration);
    if(ReadConf.fail()){
        std::cerr << "Initial configuration file '"<<initial_configuration 
                  << "' not found. Aborting.\n";
        exit(1);
    }
    std::cout << "Classic Lennard-Jones fluid\n";
    std::cout << "Monte Carlo simulation\n\n";
    std::cout << "Interatomic potential v(r) = 4 * [(1/r)^12 - (1/r)^6]\n\n";
    std::cout << "Boltzmann weight exp(- beta * sum_{i<j} v(r_ij) ), beta = 1/T" 
                 "\n\n";
    std::cout << "The program uses Lennard-Jones units\n" << std::endl;

    ReadInput >> temp;
    beta = 1.0/temp;
    std::cout << "Temperature = " << temp << std::endl;
    
    ReadInput >> npart;
    std::cout << "Number of particles = " << npart << std::endl;
    
    ReadInput >> rho;
    std::cout << "Density of particles = " << rho << std::endl;
    vol = (double)npart/rho;
    box = std::pow(vol,1.0/3.0);
    std::cout << "Volume of the simulation box = " << vol << std::endl;
    std::cout << "Edge of the simulation box = " << box << std::endl;
    
    ReadInput >> rcut;
    std::cout << "Cutoff of the interatomic potential = " << rcut 
              <<std:: endl << std::endl;
      
    //Tail corrections for potential energy and pressure
    vtail = (8.0*M_PI*rho)/(9.0*std::pow(rcut,9)) - (8.0*M_PI*rho)/
        (3.0*std::pow(rcut,3));
    ptail = (32.0*M_PI*rho)/(9.0*std::pow(rcut,9)) - (16.0*M_PI*rho)/
        (3.0*std::pow(rcut,3));
    std::cout << "Tail correction for the potential energy = " 
              << vtail << std::endl;
    std::cout << "Tail correction for the virial           = " 
              << ptail << std::endl; 
    
    ReadInput >> delta;
    
    ReadInput >> nblk;
    
    ReadInput >> nstep;
    
    std::cout << "The program perform Metropolis moves with uniform"
       " translations" << std::endl;
    std::cout << "Moves parameter = " << delta <<std:: endl;
    std::cout << "Number of blocks = " << nblk << std::endl;
    std::cout << "Number of steps in one block = " << nstep << std::endl 
              << std::endl;
    ReadInput.close();

    //set to zero blocking variables
    for (auto &it : keys){
        walker[it]=0.;
        block_average[it]=0.;
        global_average[it]=0.; 
        global_average2[it]=0.;
    }
    
    //Read initial configuration
    std::cout << "Read initial configuration from file " 
              << initial_configuration << std::endl << std::endl;
   
    double tmpx, tmpy, tmpz;
    x.resize(npart);
    y.resize(npart);
    z.resize(npart);
    for (unsigned int i=0; i<npart; ++i){
        ReadConf >> tmpx >> tmpy >> tmpz;
        x[i] = Pbc( tmpx * box );
        y[i] = Pbc( tmpy * box );
        z[i] = Pbc( tmpz * box );
    }
    ReadConf.close();
   
    bin_size = (box/2.0)/(double)nbins;
    
    std::ofstream binning("outputs/binning.dat");
    for(unsigned int i=0; i<nbins; ++i)
        binning << i*bin_size+bin_size/2. << std::endl;
    binning.close();

    histo_walker.resize(nbins);
    histo_block_average.resize(nbins);
    histo_global_average.resize(nbins);
    histo_global_average2.resize(nbins);

    //Evaluate potential energy and virial of the initial configuration
    Measure();
    //Print initial values for the potential energy and virial
    std::cout << "Initial potential energy (with tail corrections) = " 
              << walker.at("energy")/(double)npart + vtail << std::endl;
    std::cout << "Virial                   (with tail corrections) = " 
              << walker.at("virial")/(double)npart + ptail << std::endl;
    std::cout << "Pressure                 (with tail corrections) = " 
              << rho*temp+(walker.at("virial")+(double)npart*ptail)/vol 
              << std::endl << std::endl;


    blk_norm=0;
    attempted=0;
    accepted=0;
    Epot.open("outputs/epot.dat"); //,std::ios::app);
    Pres.open("outputs/pres.dat"); //,std::ios::app);
    Gerr.open("outputs/gerr.dat"); //,std::ios::app);
    Gave.open("outputs/gave.dat"); //,std::ios::app);

    if(istant){
        ist_pot.open("outputs/istant_energy.dat");
        ist_pres.open("outputs/istant_press.dat");
    }
}

double MolecularMC::Boltzmann(double xx, double yy, double zz, unsigned int ip)
{
  double ene=0.0;
  double dx, dy, dz, dr;

  for (unsigned int i=0; i<npart; ++i)
  {
    if(i != ip)
    {
// distance ip-i in pbc
      dx = Pbc(xx - x[i]);
      dy = Pbc(yy - y[i]);
      dz = Pbc(zz - z[i]);

      dr = dx*dx + dy*dy + dz*dz;
      dr = sqrt(dr);

      if(dr < rcut)
      {
        ene += 1.0/pow(dr,12) - 1.0/pow(dr,6);
      }
    }
  }

  return 4.0*ene;
}

void MolecularMC::ConfFinal(){
    std::ofstream WriteConf;
    
    std::cout << "Print final configuration to file config.final " 
              << std::endl << std::endl;
    WriteConf.open("config.final");
    for (unsigned int i=0; i<npart; ++i)
        WriteConf<<x[i]/box<<"\t"<<y[i]/box<<"\t"<<z[i]/box<<std::endl;
    WriteConf.close();
    
    rand.SaveSeed();
}

void MolecularMC::ConfXYZ(unsigned int nconf){ 
    std::ofstream WriteXYZ;
    
    WriteXYZ.open("frames/config_" + std::to_string(nconf) + ".xyz");
    WriteXYZ << npart << std::endl;
    WriteXYZ << "This is only a comment!" << std::endl;
    for (unsigned int i=0; i<npart; ++i)
        WriteXYZ << "LJ  " << Pbc(x[i]) << "   " <<  Pbc(y[i]) << "   " << 
            Pbc(z[i]) << std::endl;
    WriteXYZ.close();
}

void MolecularMC::Accumulate(){

    for(auto &it : keys)
        block_average.at(it) = block_average.at(it) + walker.at(it); 
    
    for(unsigned int i=0; i<nbins; ++i)
        histo_block_average[i]+=histo_walker[i];

    blk_norm++;
}

void MolecularMC::Reset(unsigned int iblk) //Reset block averages
{
    //TODO useful??? Don't think so...
    if(iblk==0){
        for(auto &it : keys){
            global_average.at(it) = 0.;
            global_average2.at(it) = 0.;
        }
       
        for(auto &it: histo_global_average)
            it=0.;
        for(auto &it: histo_global_average2)
            it=0.;
    }
   
    for(auto &it : keys)
        block_average.at(it) = 0.;
   
    for(auto &it: histo_block_average)
        it=0.;

    blk_norm = 0.;
    attempted = 0.;
    accepted = 0.;
}

void MolecularMC::Run(){
    unsigned int nconf = 1;
    for(unsigned int iblk=1; iblk <= nblk; ++iblk) //Simulation
    {
        Reset(iblk);   //Reset block averages
        for(unsigned int istep=1; istep <= nstep; ++istep)
        {
            Move();
            Measure();
            Accumulate(); //Update block averages
            
            if(istep%10 == 0){
                ConfXYZ(nconf);
                nconf++; 
            }
        }
        Averages(iblk);   //Print results for current block
    }
    ConfFinal(); //Write final configuration
}

void MolecularMC::Averages(unsigned int iblk) //Print results for current block
{
    
    //double r, gdir;
//    std::ofstream Gerr, Gave, Epot, Pres;
    const int wd=12;
    
    std::cout << "Block number " << iblk << std::endl;
    std::cout << "Acceptance rate " << (double)accepted/attempted << std::endl 
              << std::endl;
    
  
 
    double stima_pot, err_pot;
    stima_pot = block_average.at("energy")/blk_norm/(double)npart + vtail; 
    global_average.at("energy") += stima_pot;
    global_average2.at("energy") += stima_pot*stima_pot;
    err_pot=Error(global_average.at("energy"),
            global_average2.at("energy"),iblk);
    
    double stima_pres, err_press;
    stima_pres = rho*temp+(block_average.at("virial")/blk_norm+
            ptail*(double)npart)/vol; 
    global_average.at("virial") += stima_pres;
    global_average2.at("virial") += stima_pres*stima_pres;
    err_press=Error(global_average.at("virial"),
            global_average2.at("virial"),iblk);

    double stima_g; //, err_g;
    std::vector<double> err_g;
    for(unsigned int i=0; i<nbins; ++i){
        stima_g = (double)histo_block_average[i]/blk_norm;
        histo_global_average[i] += stima_g;
        histo_global_average2[i] += stima_g*stima_g;
        err_g.push_back(Error(histo_global_average[i],
                    histo_global_average2[i],iblk));
    }
    Gerr << std::endl;

    //Potential energy per particle
    Epot << std::setw(wd) << iblk <<std::setw(wd) << stima_pot << std::setw(wd) 
         << global_average.at("energy")/(double)iblk << std::setw(wd) << err_pot 
         << std::endl;
    //Pressure
    Pres << std::setw(wd) << iblk <<std::setw(wd) << stima_pres << std::setw(wd) 
         << global_average.at("virial")/(double)iblk << std::setw(wd) 
         << err_press << std::endl;

//g(r)
/*
    double sum=0;
    for(auto & it : histo_global_average)
        sum+=it; 
    sum=1.;
    for(auto & it : histo_global_average)
        Gave << it/sum << "\t";
    Gave << std::endl;
    for(auto & it : err_g)
        Gerr << it/sum << "\t";
    Gerr << std::endl;
*/
    double norm=0.; // *(4.*M_PI/3.);
    for(auto & it : histo_global_average)
        norm+=it; 
    for(auto & it : histo_global_average)
        Gave << it/(norm) << "\t";
    Gave << std::endl;
    for(auto & it : err_g)
        Gerr << it/(norm)<< "\t";
    Gerr << std::endl;

    std::cout << "----------------------------" << std::endl << std::endl;
}

void MolecularMC::Measure(){
    //unsigned int bin;
    double v = 0.0, w = 0.0;
    double vij, wij;
    double dx, dy, dz, dr;

    //reset the hystogram of g(r)
    for(auto &it: histo_walker)
        it=0;

    //cycle over pairs of particles
    for (unsigned int i=0; i<npart-1; ++i){
        for (unsigned int j=i+1; j<npart; ++j){

            // distance i-j in pbc
            dx = Pbc(x[i] - x[j]);
            dy = Pbc(y[i] - y[j]);
            dz = Pbc(z[i] - z[j]);

            dr = dx*dx + dy*dy + dz*dz;
            dr = sqrt(dr);
//            double r;
            //update of the histogram of g(r)
            for (unsigned int k=0; k<nbins; ++k){
                if(dr>bin_size*k and dr<=bin_size*(k+1)){
  //                  r=std::sqrt(x[j]*x[j]+y[j]*y[j]+z[j]*z[j]);
                    histo_walker[k]+=3./(std::pow(dr+bin_size,3)-std::pow(dr,3))
                        /(2.*M_PI);
                    break;
                }
            }

            if(dr < rcut){
                vij = 1.0/std::pow(dr,12) - 1.0/std::pow(dr,6);
                wij = 1.0/std::pow(dr,12) - 0.5/std::pow(dr,6);

                // contribution to energy and virial
                v += vij;
                w += wij;
            }
        }          
    } 
    walker.at("energy") = 4.0 * v;
    walker.at("virial") = 48.0 * w / 3.0;

    if (istant){
        ist_pot << walker.at("energy")/(double)npart + vtail << std::endl;
        ist_pres <<  rho*temp+(walker.at("virial")+ ptail*(double)npart)/vol
                   << std::endl;
    }
}

void MolecularMC::Move(){
    int o;
    double p, energy_old, energy_new;
    double xold, yold, zold, xnew, ynew, znew;
    
    for(unsigned int i=0; i<npart; ++i){
        //Select randomly a particle 
        o = (int)(rand.Rannyu()*npart);
        //Old
        xold = x[o];
        yold = y[o];
        zold = z[o];
        
        energy_old = Boltzmann(xold,yold,zold,o);
        
        //New
        xnew = Pbc( x[o] + delta*(rand.Rannyu() - 0.5) );
        ynew = Pbc( y[o] + delta*(rand.Rannyu() - 0.5) );
        znew = Pbc( z[o] + delta*(rand.Rannyu() - 0.5) );
        
        energy_new = Boltzmann(xnew,ynew,znew,o);
       
        //Metropolis test
        p = exp(beta*(energy_old-energy_new));
        if(p >= rand.Rannyu()){
            //Update
            x[o] = xnew;
            y[o] = ynew;
            z[o] = znew;
            
            accepted++;
        }
        attempted++;
    }
}


MolecularMC::~MolecularMC(){
    Epot.close();
    Pres.close();
    Gerr.close(); 
    Gave.close();
    if(istant){
        ist_pot.close();
        ist_pres.close();
    }
}
