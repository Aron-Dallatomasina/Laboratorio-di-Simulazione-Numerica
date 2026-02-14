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
#include <cmath>
#include <iomanip>
#include "MD_MC.h"

using namespace std;

/*il codice risolve, con simulazione, la simulaz di liq/gas/sol di un gas nobile
- in MD: sto studiando nel microcanonico (sistema isolato, che cons energia)
- in MC: sto studiando nel canonico (sistema con numero particelle, volume e T costante) [fissata non è E ma T]
	in questo caso campioneremo il peso statistico di Boltzmann col Metropolis
*/


int main() { 
	Input();										// Initialization and input handling
	int nconf=1;									// Configuration index (used for optional XYZ output)
	for(int iblk=1; iblk <= nblk; iblk++) {			// Simulation starts here - Loop over blocks
		Reset(iblk);								// Reset block statistics before starting a new block
		for(int istep=1; istep<=nstep; istep++) {	// Loop over steps within the block
			Move();									// Perform one MC (Metropolis) or MD (Verlet) move
			Measure();								// Measure instantaneous values of observables
			Accumulate();							// Accumulate observables for block averages
			if(istep%10==0){
			//ConfXYZ(nconf);						// Optional: write configuration in XYZ format. Commented to avoid "filesystem full"!
			nconf+=1;
			}
		}
		Averages(iblk);								// Compute block averages and update global statistics
	}
	ConfFinal();									// Save final configuration

	return 0;
}



// Read input parameters and initialize the system
void Input(void) {

	cout<<"Classic Lennard-Jones fluid        "<<endl;
	cout<<"MD(NVE) / MC(NVT) simulation       "<<endl<<endl;
	cout<<"Interatomic potential v(r) = 4 * [(1/r)^12 - (1/r)^6]"<<endl<<endl;
	cout<<"Boltzmann weight exp(- beta * sum_{i<j} v(r_ij) ), beta = 1/T "<<endl<<endl;
	cout<<"The program uses Lennard-Jones units "<<endl;

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
	ReadInput.open("../Input/input.in");
	ReadInput>>iNVET;							// 0=MD, 1=MC
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

	// Read simulation parameters from input file (reduced units): temperature, number of particles, density, volume, box length, cutoff radius
	ReadInput>>temp;
	beta=1.0/temp;
	cout<<"Temperature = "<<temp<<endl;
	ReadInput>>npart;
	cout<<"Number of particles = "<<npart<<endl;
	ReadInput>>rho;
	cout<<"Density of particles = "<<rho<<endl;
	vol=(double)npart/rho;
	box=pow(vol,1.0/3.0);
	cout<<"Volume of the simulation box = "<<vol<<endl;
	cout<<"Edge of the simulation box = "<<box<<endl;
	ReadInput>>rcut;
	cout<<"Cutoff of the interatomic potential = "<<rcut<<endl<<endl;

	// Read simulation parameters from input file: time step (MD) or maximum displacement (MC), number of blocks, steps per block
	ReadInput>>delta;
	ReadInput>>nblk;
	ReadInput>>nstep;
	cout<<"The program performs Metropolis moves with uniform translations"<<endl;
	cout<<"Moves parameter = "<<delta<<endl;
	cout<<"Number of blocks = "<<nblk<<endl;
	cout<<"Number of steps in one block = "<<nstep<<endl<<endl;
	ReadInput.close();

	// Define indices of the observables in the walker array
	iv=0;						// Potential energy
	it=1;						// Temperature
	ik=2;						// Kinetic energy
	ie=3;						// Total energy
	ip=4;						// Pressure (to be computed via virial theorem)
	n_props=5;					// Total number of observables

	//Read initial configuration
	ifstream ReadConf, ReadVelocity;
	cout<<"Read initial configuration"<<endl<<endl;
	if(restart) {								// Restart from previous simulation: read final configuration and velocities
		ReadConf.open("../Output/config.out");
		ReadVelocity.open("../Output/velocity.out");
		for (int i=0; i<npart; ++i) ReadVelocity>>vx[i]>>vy[i]>>vz[i];
	} else {									// Start from input configuration and generate velocities from Maxwell-Boltzmann distribution
		ReadConf.open("../Input/config.in");
		cout<<"Prepare velocities with center of mass velocity equal to zero "<<endl;
		double sumv[3]={0.0,0.0,0.0};
		for (int i=0; i<npart; ++i) {
			vx[i]=random.Gauss(0.,sqrt(temp));
			vy[i]=random.Gauss(0.,sqrt(temp));
			vz[i]=random.Gauss(0.,sqrt(temp));
			sumv[0]+=vx[i];
			sumv[1]+=vy[i];
			sumv[2]+=vz[i];
		}
		for (int idim=0; idim<3; ++idim) sumv[idim]/=(double)npart;		// Center-of-mass velocity
		double sumv2=0.0;
		for (int i=0; i<npart; ++i) {
			vx[i]=vx[i]-sumv[0];
			vy[i]=vy[i]-sumv[1];
			vz[i]=vz[i]-sumv[2];
			sumv2+=vx[i]*vx[i]+vy[i]*vy[i]+vz[i]*vz[i];
		}
		sumv2/=(double)npart;											// Mean squared velocity (after center-of-mass velocity removal)
		// Rescale velocities to match target temperature
		double fs=sqrt(3*temp/sumv2);
		cout<<"velocity scale factor: "<<fs<<endl<<endl;
		for (int i=0; i<npart; ++i) {
			vx[i]*=fs;
			vy[i]*=fs;
			vz[i]*=fs;
		}
	}

	// Read particle positions, convert to reduced units, apply periodic boundary conditions
	for (int i=0; i<npart; ++i) {												
		ReadConf>>x[i]>>y[i]>>z[i];
		x[i]=Pbc(x[i]*box);
		y[i]=Pbc(y[i]*box);
		z[i]=Pbc(z[i]*box);
	}
	ReadConf.close();

	// Initialize previous positions for Verlet integration
	for (int i=0; i<npart; ++i) {
		if(iNVET) {
			xold[i]=x[i];
			yold[i]=y[i];
			zold[i]=z[i];
		} else {
			xold[i]=Pbc(x[i]-vx[i]*delta);
			yold[i]=Pbc(y[i]-vy[i]*delta);
			zold[i]=Pbc(z[i]-vz[i]*delta);
		}
	}

	// Evaluate properties of the initial configuration
	Measure();

	//Print initial values for measured properties
	cout<<"Initial potential energy = "<<walker[iv]/(double)npart<<endl;
	cout<<"Initial temperature      = "<<walker[it]<<endl;
	cout<<"Initial kinetic energy   = "<<walker[ik]/(double)npart<<endl;
	cout<<"Initial total energy     = "<<walker[ie]/(double)npart<<endl;
	cout<<"Initial pressure         = "<<walker[ip]<<endl;

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

	ofstream Epot, Ekin, Etot, Temp, Pres;
	const int wd=12;									// Output column width
	cout<<"Block number "<<iblk<<endl;
	cout<<"Acceptance rate "<<accepted/attempted<<endl<<endl;
	Epot.open("../Output/output_epot.dat",ios::app);
	Ekin.open("../Output/output_ekin.dat",ios::app);
	Temp.open("../Output/output_temp.dat",ios::app);
	Etot.open("../Output/output_etot.dat",ios::app);
	Pres.open("../Output/output_pres.dat",ios::app);

	// Potential energy per particle
	stima_pot=blk_av[iv]/blk_norm/(double)npart;
	glob_av[iv]+=stima_pot;
	glob_av2[iv]+=stima_pot*stima_pot;
	err_pot=Error(glob_av[iv],glob_av2[iv],iblk);
	// Kinetic energy per particle
	stima_kin=blk_av[ik]/blk_norm/(double)npart;
	glob_av[ik]+=stima_kin;
	glob_av2[ik]+=stima_kin*stima_kin;
	err_kin=Error(glob_av[ik],glob_av2[ik],iblk);
	// Total energy per particle
	stima_etot=blk_av[ie]/blk_norm/(double)npart;
	glob_av[ie]+=stima_etot;
	glob_av2[ie]+=stima_etot*stima_etot;
	err_etot=Error(glob_av[ie],glob_av2[ie],iblk);
	// Temperature
	stima_temp=blk_av[it]/blk_norm;
	glob_av[it]+=stima_temp;
	glob_av2[it]+=stima_temp*stima_temp;
	err_temp=Error(glob_av[it],glob_av2[it],iblk);
	// Pressure 
	stima_pres=blk_av[ip]/blk_norm;
	glob_av[ip]+=stima_pres;
	glob_av2[ip]+=stima_pres*stima_pres;
	err_press=Error(glob_av[ip],glob_av2[ip],iblk);

	// Write block index, block estimate, progressive average and error
	Epot<<setw(wd)<<iblk<< setw(wd)<<stima_pot<<setw(wd)<<glob_av[iv]/(double)iblk<<setw(wd)<<err_pot<<endl;
	Ekin<<setw(wd)<<iblk<< setw(wd)<<stima_kin<<setw(wd)<<glob_av[ik]/(double)iblk<<setw(wd)<<err_kin<<endl;
	Etot<<setw(wd)<<iblk<< setw(wd)<<stima_etot<<setw(wd)<<glob_av[ie]/(double)iblk<<setw(wd)<<err_etot<<endl;
	Temp<<setw(wd)<<iblk<< setw(wd)<<stima_temp<<setw(wd)<<glob_av[it]/(double)iblk<<setw(wd)<<err_temp<<endl;
	Pres<<setw(wd)<<iblk<<setw(wd)<<stima_pres<<setw(wd)<<glob_av[ip]/(double)iblk<<setw(wd)<<err_press<<endl;

	cout<<"----------------------------"<<endl<<endl;
	Epot.close();
	Ekin.close();
	Etot.close();
	Temp.close();
	Pres.close();
}

// Perform one MC or MD move
void Move() {
	int o;													// Index of selected particle (MC)
	double p, energy_old, energy_new;
	double xnew, ynew, znew;

	if(iNVET) {												// Monte Carlo (NVT) move
		for(int i=0; i<npart; ++i) {
			o=(int)(random.Rannyu()*npart);					// Randomly select a particle (0 <= o <= npart-1)
			energy_old=Boltzmann(x[o],y[o],z[o],o);			// Energy before trial move
			// Propose random displacement within [-delta/2, +delta/2]
			x[o]=Pbc(x[o]+delta*(random.Rannyu()-0.5));
			y[o]=Pbc(y[o]+delta*(random.Rannyu()-0.5));
			z[o]=Pbc(z[o]+delta*(random.Rannyu()-0.5));
			// Evaluate
			energy_new=Boltzmann(x[o],y[o],z[o],o);			// Energy after trial move
			p=exp(beta*(energy_old-energy_new));			// Metropolis acceptance probability
			if(p>=random.Rannyu()) {						// Accept move
				xold[o]=x[o];
				yold[o]=y[o];
				zold[o]=z[o];
				accepted=accepted+1.0;
			} else {										// Reject move
				x[o]=xold[o];
				y[o]=yold[o];
				z[o]=zold[o];
			}
			attempted=attempted+1.0;
		}
	} else {												// Molecular Dynamics (NVE) move
		double fx[m_part], fy[m_part], fz[m_part];
		for(int i=0; i<npart; ++i){							// Compute forces on all particles
			fx[i]=Force(i,0);
			fy[i]=Force(i,1);
			fz[i]=Force(i,2);
		}
		for(int i=0; i<npart; ++i){							//Verlet integration scheme
			xnew=Pbc(2.0*x[i]-xold[i]+fx[i]*pow(delta,2));
			ynew=Pbc(2.0*y[i]-yold[i]+fy[i]*pow(delta,2));
			znew=Pbc(2.0*z[i]-zold[i]+fz[i]*pow(delta,2));
			vx[i]=Pbc(xnew-xold[i])/(2.0*delta);
			vy[i]=Pbc(ynew-yold[i])/(2.0*delta);
			vz[i]=Pbc(znew-zold[i])/(2.0*delta);
			// Update positions
			xold[i]=x[i];
			yold[i]=y[i];
			zold[i]=z[i];
			x[i]=xnew;
			y[i]=ynew;
			z[i]=znew;
			accepted=accepted+1.0;
			attempted=attempted+1.0;
		}
	}
	return;
}

// Save final configuration and velocities
void ConfFinal(void) {
	ofstream WriteConf, WriteVelocity, WriteSeed;
	cout<<"Print final configuration to file config.out"<<endl<<endl;
	WriteConf.open("../Output/config.out");
	WriteVelocity.open("../Output/velocity.out");
	for (int i=0; i<npart; ++i) {
		WriteConf<<x[i]/box<<"   "<<y[i]/box<<"   "<<z[i]/box<<endl;	// Reduced positions
		WriteVelocity<<vx[i]<<"   "<<vy[i]<<"   "<<vz[i]<<endl;			// Velocities
	}
	WriteConf.close();
	WriteVelocity.close();
	random.SaveSeed();							// Save RNG state for possible restart
}

// Write configuration in XYZ format
void ConfXYZ(int nconf) {
	ofstream WriteXYZ;
	WriteXYZ.open("../Output/frames/config_"+to_string(nconf)+".xyz");
	WriteXYZ<<npart<<endl;
	WriteXYZ<<"This is only a comment!"<<endl;
	for (int i=0; i<npart; ++i) WriteXYZ<<"LJ  "<<Pbc(x[i])<<"   "<<Pbc(y[i])<<"   "<<Pbc(z[i])<<endl;
	WriteXYZ.close();
}

// Measure instantaneous physical observables
void Measure() {
	double v=0.0, kin=0.0, w=0.0;					// Potential energy, kinetic energy, virial
	double vij;
	double dx,dy,dz,dr;
	// Loop over all distinct particle pairs (i<j)
	for (int i=0; i<npart-1; ++i) {
		for (int j=i+1; j<npart; ++j) {
		// distance i-j in pbc
			dx=Pbc(x[i]-x[j]);
			dy=Pbc(y[i]-y[j]);
			dz=Pbc(z[i]-z[j]);
			dr=dx*dx+dy*dy+dz*dz;
			dr=sqrt(dr);
		// Apply cutoff radius
			if(dr<rcut) {
				vij=1.0/pow(dr,12)-1.0/pow(dr,6);	// Lennard-Jones pair potential (reduced units)
				v+=vij;								// Accumulate pair contribution to potential energy
				w+=1.0/pow(dr,12)-0.5/pow(dr,6);	// Accumulate pair contribution to virial

			}
		}
	}
	// Compute kinetic energy
	for (int i=0; i<npart; ++i) kin+=0.5*(vx[i]*vx[i]+vy[i]*vy[i]+vz[i]*vz[i]);

	walker[iv]=4.0*v;								// Potential energy
	walker[ik]=kin;									// Kinetic energy
	walker[it]=(2.0/3.0)*kin/(double)npart;			// Temperature
	walker[ie]=4.0*v+kin;							// Total energy;
	walker[ip]=rho*walker[it]+16.0/vol*w;			// Pressure

	return;
}

// Compute single-particle potential energy (MC) - ip=particle index
double Boltzmann(double xx, double yy, double zz, int ip) {
	double ene=0.0;
	double dx,dy,dz,dr;
	// Sum interactions between particle ip and all other particles
	for (int i=0; i<npart; ++i) {
		if(i!=ip) {
			// distance ip-i in pbc
			dx=Pbc(xx-x[i]);
			dy=Pbc(yy-y[i]);
			dz=Pbc(zz- z[i]);
			dr=dx*dx+dy*dy+dz*dz;
			dr=sqrt(dr);
			if(dr<rcut) ene+=1.0/pow(dr,12)-1.0/pow(dr,6);		// Apply cutoff radius, use Lennard-Jones pair potential
		}
	}
	return 4.0*ene;						// Return potential energy contribution (reduced units)
}

// Apply periodic boundary conditions (side L=box)
double Pbc(double r)  {
	return r-box*rint(r/box);
}

// Statistical uncertainty (data blocking)
double Error(double sum, double sum2, int iblk) {
	return sqrt(fabs(sum2/(double)iblk-pow(sum/(double)iblk,2))/(double)iblk);
}

// Compute force component acting on a particle - ip=particle index, idir=cartesian coordinate index
double Force(int ip, int idir){
	double f=0.0;
	double dvec[3],dr;
	// Sum contributions from all other particles
	for (int i=0; i<npart; ++i){
		if(i!=ip){
			// distance ip-i in pbc
			dvec[0]=Pbc(x[ip]-x[i]);
			dvec[1]=Pbc(y[ip]-y[i]);
			dvec[2]=Pbc(z[ip]-z[i]);
			dr=dvec[0]*dvec[0]+dvec[1]*dvec[1]+dvec[2]*dvec[2];
			dr=sqrt(dr);
			// Apply cutoff radius
			if(dr<rcut) f+=dvec[idir]*(48.0/pow(dr,14)-24.0/pow(dr,8));		// Lennard-Jones force component, -Grad_ip V(r)
		}
	}
	return f;
}