#include "molecular_dynamics.h"
#include <iostream>
#include "functions.h"          //to initialize random


MolecularDynamics::MolecularDynamics(std::string simParameters, 
                                     std::string configFile)
{
    Epot.open("results/output_epot.dat");
    Ekin.open("results/output_ekin.dat");
    Temp.open("results/output_temp.dat");
    Etot.open("results/output_etot.dat");
    Press.open("results/output_press.dat");
    
    initialize_seed(rand, "Primes", "seed.in"); //seed of random initialization

    std::ifstream ReadInput(simParameters);
    if(ReadInput.fail()){
        std::cerr << "Unable to open setup file\n";
        exit(1);
    }

    std::cout << "Classic Lennard-Jones fluid\n"; 
    std::cout << "Molecular dynamics simulation in NVE ensemble  \n\n";
    std::cout << "Interatomic potential v(r) = 4 * [(1/r)^12 - (1/r)^6]\n\n";
    std::cout << "The program uses Lennard-Jones units \n";
  
    ReadInput >> temp;
    ReadInput >> npart;
    std::cout << "Number of particles = " << npart << std::endl;
    //resize vectors
    x.resize(npart);
    y.resize(npart);
    z.resize(npart);
    xold.resize(npart);
    yold.resize(npart);
    zold.resize(npart);
    vx.resize(npart);
    vy.resize(npart);
    vz.resize(npart);
    fx.resize(npart);
    fy.resize(npart);
    fz.resize(npart);

    ReadInput >> rho;
    std::cout << "Density of particles = " << rho << std::endl;
    vol = (double)npart/rho;
    std::cout << "Volume of the simulation box = " << vol << std::endl;
    box = pow(vol,1.0/3.0);
    std::cout << "Edge of the simulation box = " << box << std::endl;
  
    ReadInput >> rcut;
    ReadInput >> delta;
    ReadInput >> nstep;
    ReadInput >> iprint;
  
    std::cout
    <<  "The program integrates Newton equations with the Verlet method\n"
    << "Time step = " << delta << std::endl
    << "Number of steps = " << nstep << std::endl;

    ReadInput >> measure_time_interval;
    std::cout << "Measures performed every " << measure_time_interval << 
        "time steps.\n";
    unsigned int n_blocks;
    ReadInput >> n_blocks;
    std::cout <<"Statistical error computed with "<< n_blocks <<" blocks.\n\n";
    est_pot.resize(n_blocks);
    est_kin.resize(n_blocks);
    est_etot.resize(n_blocks);
    est_temp.resize(n_blocks);
    est_press.resize(n_blocks);
    std::fill(est_pot.begin(), est_pot.end(), 0.);
    std::fill(est_kin.begin(), est_kin.end(), 0.);
    std::fill(est_etot.begin(), est_etot.end(), 0.);
    std::fill(est_temp.begin(), est_temp.end(), 0.);
    std::fill(est_press.begin(), est_press.end(), 0.);
    

    ReadInput.close();
    //compute block size and setting block index
    block_size=(nstep/measure_time_interval)/n_blocks;
    iblock=0;
    imeasure=0;


    //Prepare array for measurements
    //Read initial configuration
    std::cout << "Read initial configuration from file "+configFile <<
        std::endl<<std::endl; 
    
    std::ifstream ReadConf(configFile);
    if(ReadConf.fail()){
        std::cerr << "Error opening config file\n";
        exit(1);
    }
    
    for (unsigned int i=0; i<npart; ++i){
        ReadConf >> x[i] >> y[i] >> z[i];
        x[i] = x[i] * box;
        y[i] = y[i] * box;
        z[i] = z[i] * box;
    }
    ReadConf.close();
  
    //Prepare initial velocities
    std::cout << "Prepare random velocities with center of mass velocity equal"
        "to zero \n\n"; 
    double sumv[3] = {0.0, 0.0, 0.0};
    for (unsigned int i=0; i<npart; ++i){
        vx[i] = rand.Rannyu() - 0.5;
        vy[i] = rand.Rannyu() - 0.5;
        vz[i] = rand.Rannyu() - 0.5;
        
        sumv[0] += vx[i];
        sumv[1] += vy[i];
        sumv[2] += vz[i];
    }
    for (int idim=0; idim<3; ++idim) sumv[idim] /= (double)npart;
    double sumv2 = 0.0, fs;
    for (unsigned int i=0; i<npart; ++i){
        vx[i] = vx[i] - sumv[0];
        vy[i] = vy[i] - sumv[1];
        vz[i] = vz[i] - sumv[2];
    
        sumv2 += vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i];
    }
    sumv2 /= (double)npart;
    
    fs = sqrt(3 * temp / sumv2);   // fs = velocity scale factor 
    for (unsigned int i=0; i<npart; ++i){
        vx[i] *= fs;
        vy[i] *= fs;
        vz[i] *= fs;
        xold[i] = x[i] - vx[i] * delta;
        yold[i] = y[i] - vy[i] * delta;
        zold[i] = z[i] - vz[i] * delta;
    }


}

MolecularDynamics::MolecularDynamics(std::string simParameters, 
                                     std::string configFile,
                                     std::string oldConfigFile)
{
    Epot.open("results/output_epot.dat");
    Ekin.open("results/output_ekin.dat");
    Temp.open("results/output_temp.dat");
    Etot.open("results/output_etot.dat");
    Press.open("results/output_press.dat");
 
    initialize_seed(rand, "Primes", "seed.in"); //seed of random initialization
    
    std::ifstream ReadInput(simParameters);
    if(ReadInput.fail()){
        std::cerr << "Unable to open setup file\n";
        exit(1);
    }

    std::cout << "Classic Lennard-Jones fluid\n"; 
    std::cout << "Molecular dynamics simulation in NVE ensemble  \n\n";
    std::cout << "Interatomic potential v(r) = 4 * [(1/r)^12 - (1/r)^6]\n\n";
    std::cout << "The program uses Lennard-Jones units \n";
  
    ReadInput >> temp;
    ReadInput >> npart;
    std::cout << "Number of particles = " << npart << std::endl;
    //resize vectors
    x.resize(npart);
    y.resize(npart);
    z.resize(npart);
    xold.resize(npart);
    yold.resize(npart);
    zold.resize(npart);
    vx.resize(npart);
    vy.resize(npart);
    vz.resize(npart);
    fx.resize(npart);
    fy.resize(npart);
    fz.resize(npart);

    ReadInput >> rho;
    std::cout << "Density of particles = " << rho << std::endl;
    vol = (double)npart/rho;
    std::cout << "Volume of the simulation box = " << vol << std::endl;
    box = pow(vol,1.0/3.0);
    std::cout << "Edge of the simulation box = " << box << std::endl;
  
    ReadInput >> rcut;
    ReadInput >> delta;
    ReadInput >> nstep;
    ReadInput >> iprint;
  
    std::cout
    <<  "The program integrates Newton equations with the Verlet method\n"
    << "Time step = " << delta << std::endl
    << "Number of steps = " << nstep << std::endl;
    ReadInput >> measure_time_interval;
    std::cout << "Measures performed every " << measure_time_interval << 
        "time steps.\n";
    unsigned int n_blocks;
    ReadInput >> n_blocks;
    std::cout <<"Statistical error computed with "<< n_blocks <<" blocks.\n\n";
    //resize and fill estimators
    est_pot.resize(n_blocks);
    est_kin.resize(n_blocks);
    est_etot.resize(n_blocks);
    est_temp.resize(n_blocks);
    est_press.resize(n_blocks);
    std::fill(est_pot.begin(), est_pot.end(), 0.);
    std::fill(est_kin.begin(), est_kin.end(), 0.);
    std::fill(est_etot.begin(), est_etot.end(), 0.);
    std::fill(est_press.begin(), est_press.end(), 0.);
    
    ReadInput.close();
    
    //compute block size and setting block index
    block_size=(nstep/measure_time_interval)/n_blocks;
    iblock=0;
    imeasure=0;

    //Prepare array for measurements
    //Read initial configuration
    std::cout << "Read configuration from file "+configFile<<
        std::endl<<std::endl; 
    std::ifstream ReadConf(configFile);
    if(ReadConf.fail()){
        std::cerr << "Error opening config file\n";
        exit(1);
    }
    for (unsigned int i=0; i<npart; ++i){
        ReadConf >> x[i] >> y[i] >> z[i];
        x[i] = x[i] * box;
        y[i] = y[i] * box;
        z[i] = z[i] * box;
    }
    ReadConf.close();
    
    //read previous configuration file
    std::cout << "Read old configuration from file "+oldConfigFile<<
        std::endl<<std::endl; 
    ReadConf.open(oldConfigFile);
    if(ReadConf.fail()){
        std::cerr << "Error opening config file\n";
        exit(1);
    }
    
    for (unsigned int i=0; i<npart; ++i){
        ReadConf >> xold[i] >> yold[i] >> zold[i];
        xold[i] = xold[i] * box;
        yold[i] = yold[i] * box;
        zold[i] = zold[i] * box;
    }
    ReadConf.close();

    //one step of verlet
    Force();

    //will be new coordnate
    std::vector<double> xtmp(npart), ytmp(npart), ztmp(npart);
    for(unsigned int i=0; i<npart; ++i){ //Verlet integration scheme
        xtmp[i] = Pbc( 2.0 * x[i] - xold[i] + fx[i] * pow(delta,2) );
        ytmp[i] = Pbc( 2.0 * y[i] - yold[i] + fy[i] * pow(delta,2) );
        ztmp[i] = Pbc( 2.0 * z[i] - zold[i] + fz[i] * pow(delta,2) );
  
        vx[i] = Pbc(xtmp[i] - xold[i])/(2.0 * delta);
        vy[i] = Pbc(ytmp[i] - yold[i])/(2.0 * delta);
        vz[i] = Pbc(ztmp[i] - zold[i])/(2.0 * delta);
    }

    //compute velocities
    double t=0.;
    for (unsigned int i=0; i<npart; ++i)
        t+= (vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i]);
    stima_temp= t/(npart*3.); 
     
    double scaling_factor=std::sqrt(temp/stima_temp);
    std::cout << "scaling: " << scaling_factor << std::endl;
 
    for(unsigned int i=0; i<npart; ++i){
        vx[i]=scaling_factor*vx[i];
        vy[i]=scaling_factor*vy[i];
        vz[i]=scaling_factor*vz[i];
    }
    for(unsigned int i=0; i<npart; ++i){
        xold[i]=xtmp[i]-2*delta*vx[i];
        yold[i]=ytmp[i]-2*delta*vy[i];
        zold[i]=ztmp[i]-2*delta*vz[i];
    }
}


MolecularDynamics::~MolecularDynamics()
{
    Epot.close();
    Ekin.close();
    Temp.close();
    Etot.close();
    Press.close();
}


void MolecularDynamics::ConfFinal(std::string filename) const
{ 
    std::ofstream WriteConf;
  
    std::cout << "Print configuration to "+ filename + "\n";
    WriteConf.open(filename);
  
    for (unsigned int i=0; i<npart; ++i)
        WriteConf << x[i]/box << "   " <<  y[i]/box << "   " << z[i]/box 
            << std::endl;
  
    WriteConf.close();
}

void MolecularDynamics::ConfXYZ(unsigned int nconf) const
{ 
    std::ofstream WriteXYZ;
  
    WriteXYZ.open("frames/config_" + std::to_string(nconf) + ".xyz");
    WriteXYZ << npart << std::endl;
    WriteXYZ << "This is only a comment!" << std::endl;
    for (unsigned int i=0; i<npart; ++i)
        WriteXYZ << "LJ  " << Pbc(x[i]) << "   " <<  Pbc(y[i]) << "   " << 
            Pbc(z[i]) << std::endl;
   
    WriteXYZ.close();
}

void MolecularDynamics::Move()
{ //Move particles with Verlet algorithm
    double xnew, ynew, znew;
    //compute forces
    Force();
    for(unsigned int i=0; i<npart; ++i){ //Verlet integration scheme
        xnew = Pbc( 2.0 * x[i] - xold[i] + fx[i] * pow(delta,2) );
        ynew = Pbc( 2.0 * y[i] - yold[i] + fy[i] * pow(delta,2) );
        znew = Pbc( 2.0 * z[i] - zold[i] + fz[i] * pow(delta,2) );
        
        vx[i] = Pbc(xnew - xold[i])/(2.0 * delta);
        vy[i] = Pbc(ynew - yold[i])/(2.0 * delta);
        vz[i] = Pbc(znew - zold[i])/(2.0 * delta);
     
        xold[i] = x[i];
        yold[i] = y[i];
        zold[i] = z[i];
        
        x[i] = xnew;
        y[i] = ynew;
        z[i] = znew;
    }
}

void MolecularDynamics::Force()
{ //Compute forces as -Grad_ip V(r)
    double d1, d2,d3, dr;
    double moltiplier;
    std::fill(fx.begin(), fx.end(), 0.);
    std::fill(fy.begin(), fy.end(), 0.);
    std::fill(fz.begin(), fz.end(), 0.);
    for(unsigned int j=0; j<npart; ++j){
        for (unsigned int i=0; i<npart; ++i){
            if(i != j){
                d1 = Pbc( x[j] - x[i] );  // distance ip-i in pbc
                d2 = Pbc( y[j] - y[i] );
                d3 = Pbc( z[j] - z[i] );

                dr = d1*d1 + d2*d2 + d3*d3;
                dr = sqrt(dr);
                if(dr < rcut){
                    moltiplier=(48.0/pow(dr,14) - 24.0/pow(dr,8));
                    fx[j] += d1*moltiplier;
                    fy[j] += d2*moltiplier;
                    fz[j] += d3*moltiplier;
                }
            }
        }
    }
}


void MolecularDynamics::Measure()
{ //Properties measurement
    double v, t;
    double dx, dy, dz, dr;
    
    v = 0.0; //reset observables
    t = 0.0;
   stima_press=0.;
//cycle over pairs of particles
    for (unsigned int i=0; i<npart-1; ++i){
        for (unsigned int j=i+1; j<npart; ++j){
            dx = Pbc( x[i] - x[j] );
            dy = Pbc( y[i] - y[j] );
            dz = Pbc( z[i] - z[j] );
  
            dr = dx*dx + dy*dy + dz*dz;
            dr = sqrt(dr);
  
            if(dr < rcut){
                v+= 4.0/pow(dr,12) - 4.0/pow(dr,6);
                stima_press+=1./pow(dr,12)-0.5/pow(dr,6);
            }
         }          
    }

//Kinetic energy
    for (unsigned int i=0; i<npart; ++i)
          t += 0.5 * (vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i]);
    
    stima_pot = v/(double)npart; //Potential energy
    stima_kin = t/(double)npart; //Kinetic energy
    stima_temp = (2.0 / 3.0) * t/(double)npart; //Temperature
    stima_etot = (t+v)/(double)npart; //Total enery
    stima_press=16.*stima_press/(npart*vol);
    stima_press+=stima_temp*rho;

    Epot << stima_pot  << std::endl;
    Ekin << stima_kin  << std::endl;
    Temp << stima_temp << std::endl;
    Etot << stima_etot << std::endl;
    Press<< stima_press << std::endl;

    imeasure++;
    iblock=imeasure/block_size;
    est_pot[iblock]+=stima_pot;    
    est_kin[iblock]+=stima_kin;    
    est_temp[iblock]+=stima_temp;    
    est_etot[iblock]+=stima_etot;    
    est_press[iblock]+=stima_press;    
}

void MolecularDynamics::PrintBlocking()
{
    for(auto &it: est_pot)
        it/=block_size;

    for(auto &it: est_kin)
        it/=block_size;
    
    for(auto &it: est_temp)
        it/=block_size;

    for(auto &it: est_etot)
        it/=block_size;
    
    for(auto &it: est_press)
        it/=block_size;
    
    //blocking_error is from "function.h"
    std::vector<double> pot_err= blocking_error(est_pot);
    std::vector<double> kin_err= blocking_error(est_kin);
    std::vector<double> temp_err= blocking_error(est_temp);
    std::vector<double> etot_err= blocking_error(est_etot);
    std::vector<double> press_err= blocking_error(est_press);

    std::ofstream out("results/ave_epot.out");
    for(unsigned int i=0; i<est_pot.size(); ++i)
        out << est_pot[i] << "\t" << pot_err[i] << std::endl;
    out.close();

    out.open("results/ave_ekin.out");
    for(unsigned int i=0; i<est_kin.size(); ++i)
        out << est_kin[i] << "\t" << kin_err[i] << std::endl;
    out.close();

    out.open("results/ave_etot.out");
    for(unsigned int i=0; i<est_etot.size(); ++i)
        out << est_etot[i] << "\t" << etot_err[i] << std::endl;
    out.close();

    out.open("results/ave_temp.out");
    for(unsigned int i=0; i<est_temp.size(); ++i)
        out << est_temp[i] << "\t" << temp_err[i] << std::endl;
    out.close();
    
    out.open("results/ave_press.out");
    for(unsigned int i=0; i<est_press.size(); ++i)
        out << est_press[i] << "\t" << press_err[i] << std::endl;
    out.close();
    
}

void MolecularDynamics::RunSimulation()
{
    unsigned int nconf=1;
    std::ofstream out("frames/config_1.xyz"); //check if frames folder exists
    if(out.fail()){
        std::cout << "\n\nERROR: Unable to write on folder 'frames', create it"
            "and run again!\n\n";
        exit(2);
    }
    out.close();
    //modified in order to print one xyz file less: the corresponding 
    //configuration will be printed in the last lines to "config.old" and
    //"config.final".
    for (unsigned int istep=1; istep<nstep; ++istep){
        Move();
        if(istep % iprint ==0) 
            std::cout << "Number of time-steps: " << istep << std::endl;
        if(istep % measure_time_interval == 0){

            Measure();
            ConfXYZ(nconf);
            nconf++;
        }
    }
    std::cout << "Number of time-steps: " << nstep << std::endl << std::endl;
    ConfFinal("config.old");
    Move();
    Measure();
    ConfFinal("config.final");
    PrintBlocking();
}

