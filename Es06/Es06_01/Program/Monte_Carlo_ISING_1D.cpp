/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#include <iostream>
#include <fstream>
#include <ostream>
#include <sstream>
#include <cmath>
#include <iomanip>
#include <string>
#include "Monte_Carlo_ISING_1D.h"

using namespace std;

int main() { 
	Input();										// Initialization and input handling
	for(int iblk=1; iblk<=nblk; iblk++) {			// Simulation starts here - Loop over blocks
		Reset(iblk);								// Reset block statistics before starting a new block
		for(int istep=1; istep<=nstep; istep++) {	// Loop over steps within the block
			Move();									// Perform one move (Metropolis or Gibbs)
			Measure();								// Measure instantaneous values of observables
			Accumulate();							// Accumulate observables for block averages
		}
		Averages(iblk);								// Compute block averages and update global statistics
	}
	ConfFinal();									// Save final configuration

	return 0;
}


// Read input parameters and initialize the system
void Input(void) {

	cout<<"Classic 1D Ising model             "<<endl;
	cout<<"Monte Carlo simulation             "<<endl<<endl;
	cout<<"Nearest neighbour interaction      "<<endl<<endl;
	cout<<"Boltzmann weight exp(- beta * H ), beta = 1/T "<<endl<<endl;
	cout<<"The program uses k_B=1 and mu_B=1 units "<<endl;


	// Read the prime addend (c) used to select the random stream
	ifstream inPrime("../Input/Primes");
	if (!inPrime.is_open()) {
		cerr<<"Input error: Unable to open Primes"<<endl;
		exit(1);
	}
	int p1, p2;
	inPrime>>p1>>p2;
	inPrime.close();

	//Read input informations
	ifstream ReadInput;
	ReadInput.open("../Input/input.dat");
	ReadInput>>metro;							// 0=Gibbs, 1=Metropolis
	ReadInput>>restart;							// 0=restart the simulation, 1=continue a previous simulation

	// Read the initial seed and set the generator state
	ifstream inSeed;
	if(restart)	inSeed.open("../Output/seed.out");
	else		inSeed.open("../Input/seed.in");
	if (!inSeed.is_open()) {
		cerr<<"Input error: Unable to open seed file"<<endl;
		exit(1);
	}
	inSeed>>seed[0]>>seed[1]>>seed[2]>>seed[3];
	random.SetRandom(seed,p1,p2);
	inSeed.close();

	// Read simulation parameters from input file (reduced units): temperature, number of spins, exchange interaction, external magnetic field
	ReadInput>>temp;
	beta=1.0/temp;
	cout<<"Temperature = "<<temp<<endl;
	ReadInput>>nspin;
	cout<<"Number of spins = "<<nspin<<endl;
	ReadInput>>J;
	cout<<"Exchange interaction = "<<J<<endl;
	ReadInput>>h;
	cout<<"External field = "<<h<<endl<<endl;

	// Read simulation parameters from input file: number of blocks, steps per block
	ReadInput>>nblk;
	ReadInput>>nstep;
	if(metro==1)	cout<<"The program perform Metropolis moves"<<endl;
	else			cout<<"The program perform Gibbs moves"<<endl;
	cout<<"Number of blocks = "<<nblk<<endl;
	cout<<"Number of steps in one block = "<<nstep<<endl<<endl;
	ReadInput.close();

	// Define indices of the observables in the walker array
	iu=0;						// Energy
	ic=1;						// Heat capacity
	im=2;						// Magnetization
	ix=3;						// Magnetic susceptibility
	n_props=4;					// Total number of observables

	//Read initial configuration
	if(restart) {								// Restart from previous simulation: read final spins configuration
		ifstream ReadConf;
		ReadConf.open("../Output/config.out");
		for (int i=0; i<nspin; ++i) ReadConf>>s[i];
		ReadConf.close();
	} else {									// Generate starting spins configuration (randomly)
		for (int i=0; i<nspin; ++i)	{
			s[i]=(random.Rannyu()>=0.5)? 1:-1;
		}
	}

	// Evaluate properties of the initial configuration
	Measure();

	return;
}

// Reset block statistics
void Reset(int iblk) {
	// At the first block, initialize global averages
	if(iblk==1) {
		for(int i=0;i<n_props;++i){
			glob_av[i]=0;
			glob_av2[i]=0;
		}
	}
	for(int i=0;i<n_props;++i) blk_av[i]=0;	// Reset block averages
	blk_norm=0;								// Reset block normalization counter
	attempted=0;							// Reset number of attempted moves
	accepted=0;								// Reset number of accepted moves
}

// Accumulate observables within a block
void Accumulate(void) {
	for(int i=0; i<n_props; ++i) blk_av[i]=blk_av[i]+walker[i];
	blk_norm=blk_norm+1.0;			// Increase block normalization counter
}

// Compute block averages and statistical uncertainties
void Averages(int iblk) {

	ofstream Ene, Heat, Mag, Chi;
	const int wd=12;									// Output column width
	cout<<"Block number "<<iblk<<endl;
	if (metro==1) cout<<"Acceptance rate "<<accepted/attempted<<endl<<endl;

	ostringstream hs, Ts;
	hs<<fixed<<setprecision(2)<<h;
	Ts<<fixed<<setprecision(2)<<temp;
	// Select h directory
	string h_dir;
	cout<<h<<endl;
	if (fabs(h-0.00)<1e-6)		h_dir="../Output/h000";
	else if (fabs(h-0.02)<1e-6)	h_dir="../Output/h002";
	else						h_dir="../Output/hothers";
	// Select algorithm subdirectory
	string algo_dir=(metro==1)? "/Metropolis":"/Gibbs";
	string full_dir=h_dir+algo_dir;
	Ene.open(full_dir+"/ene_T"+hs.str()+".dat",ios::app);
	Heat.open(full_dir+"/heat_T"+hs.str()+".dat",ios::app);
	Mag.open(full_dir+"/mag_T"+hs.str()+".dat",ios::app);
	Chi.open(full_dir+"/chi_T"+hs.str()+".dat",ios::app);
	
	// Energy per particle
	stima_u=blk_av[iu]/blk_norm/(double)nspin;
	glob_av[iu]+=stima_u;
	glob_av2[iu]+=stima_u*stima_u;
	err_u=Error(glob_av[iu],glob_av2[iu],iblk);
	// Heat capacity per particle (specific heat)
	stima_c=beta*beta * ( (blk_av[ic]/blk_norm) - pow(blk_av[iu]/blk_norm,2) ) / (double)nspin;
	glob_av[ic]+=stima_c;
	glob_av2[ic]+=stima_c*stima_c;
	err_c=Error(glob_av[ic],glob_av2[ic],iblk);
	// Total magnetization per particle
	stima_m=blk_av[im]/blk_norm/(double)nspin;
	glob_av[im]+=stima_m;
	glob_av2[im]+=stima_m*stima_m;
	err_m=Error(glob_av[im],glob_av2[im],iblk);
	// Magnetic susceptibility
	stima_x=beta* ( (blk_av[ix]/blk_norm) - pow(blk_av[im]/blk_norm,2) );
	glob_av[ix]+=stima_x;
	glob_av2[ix]+=stima_x*stima_x;
	err_x=Error(glob_av[ix],glob_av2[ix],iblk);


	// Write block index, block estimate, progressive average and error
	Ene<<setw(wd)<<iblk<< setw(wd)<<stima_u<<setw(wd)<<glob_av[iu]/(double)iblk<<setw(wd)<<err_u<<endl;
	Heat<<setw(wd)<<iblk<< setw(wd)<<stima_c<<setw(wd)<<glob_av[ic]/(double)iblk<<setw(wd)<<err_c<<endl;
	Mag<<setw(wd)<<iblk<< setw(wd)<<stima_m<<setw(wd)<<glob_av[im]/(double)iblk<<setw(wd)<<err_m<<endl;
	Chi<<setw(wd)<<iblk<< setw(wd)<<stima_x<<setw(wd)<<glob_av[ix]/(double)iblk<<setw(wd)<<err_x<<endl;

	cout<<"----------------------------"<<endl<<endl;
	Ene.close();
	Heat.close();
	Mag.close();
	Chi.close();
}

void Move() {
	int o;													// Index of selected particle
	double p, energy_old, energy_new, sm;
	double energy_up, energy_down;
	for(int i=0; i<nspin; ++i) {
		o=(int)(random.Rannyu()*nspin);						// Randomly select a spin (0 <= o <= nspin-1)
		if(metro==1) {								// Metropolis
			energy_old=Boltzmann(s[o],o);					// Energy before trial move
			// Propose random spin flip
			sm=-1.0*s[o];
			// Evaluate
			energy_new=Boltzmann(sm,o);						// Energy after trial spin flip
			p=exp(beta*(energy_old-energy_new));			// Metropolis acceptance probability
			if(p>=random.Rannyu()) {						// Accept move
				s[o]=sm;
				accepted=accepted+1.0;
			}												// Reject move: don't do anything
			attempted=attempted+1.0;
		} else {									// Gibbs sampling
			energy_up=Boltzmann(1,o);
			energy_down=Boltzmann(-1,o);
			p=exp(-beta*energy_up)/(exp(-beta*energy_up)+exp(-beta*energy_down));
			if(random.Rannyu()<p)	s[o]=+1;
			else					s[o]=-1;
			attempted+=1.0;
			accepted+=1.0;									// Always accepted
		}
	}
}

// Save final spin configuration
void ConfFinal(void) {
	ofstream WriteConf;
	WriteConf.open("../Output/config.out");
	cout<<"Print final configuration to file config.out"<<endl<<endl;
	for (int i=0; i<nspin; ++i) WriteConf<<s[i]<<endl;
	WriteConf.close();
	random.SaveSeed();							// Save RNG state for possible restart
}

void Measure() {
	double u=0.0, m=0.0;
	// Cycle over spins
	for (int i=0; i<nspin; ++i) {
		u+=-J*s[i]*s[Pbc(i+1)]-0.5*h*(s[i]+s[Pbc(i+1)]);
		m+=s[i];
	}
	walker[iu]=u;		// Istantaneous energy
	walker[ic]=u*u;		// Needed for heat capacity
	walker[im]=m;		// Istantaneous magnetization
	walker[ix]=m*m;		// Needed for susceptibility
}

// Compute energy contribution of a trial spin value - ip=spin index
double Boltzmann(int sm, int ip) {
	double ene=-J*sm*(s[Pbc(ip-1)]+s[Pbc(ip+1)])-h*sm;
	return ene;
}

// Apply periodic boundary conditions
int Pbc(int i)  {
	if(i>=nspin) i=i-nspin;
	else if(i<0) i=i+nspin;
	return i;
}

// Statistical uncertainty (data blocking)
double Error(double sum, double sum2, int iblk) {
	if(iblk==1) return 0.0;
	else return sqrt((sum2/(double)iblk-pow(sum/(double)iblk,2))/(double)(iblk-1));
}