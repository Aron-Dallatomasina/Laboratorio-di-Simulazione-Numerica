/********************************************************************************************************************************
 * Numerical Simulation Laboratory - Exercises 08.1 and 08.2
 * Variational Monte Carlo and Simulated Annealing
 *
 * The program implements a Variational Monte Carlo (VMC) algorithm to estimate the ground-state energy of a one-dimensional quantum system with potential V(x) = x^4 - (5/2)x^2.
 *
 * Exercise 08.1: A trial wavefunction \Psi_{\mu,\sigma}(x), given by the sum of two Gaussians centered at ±\mu with width \sigma, is sampled using the Metropolis algorithm.
 * The expectation value of the Hamiltonian is computed via data blocking.
 *
 * Exercise 08.2: The variational parameters (\mu,\sigma) are optimized using a Simulated Annealing algorithm, where the estimated variational energy plays the role of the cost function.
 ********************************************************************************************************************************/



// C++ Standard Library
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <iomanip>
#include <cmath>
#include <string>
#include <sstream>
// C++ Standard Template Library
#include <vector>
#include <algorithm>
#include <numeric>
// RNG
#include "random.h"

constexpr double pi=3.14159265358979323846;
using namespace std;

// Metropolis move
void MetropolisMove(int& accepted, int& attempted, double& x, double& y, double delta, double mu, double sigma, Random& random);
// Return the minimum between two real numbers a and b
double Min(double a, double b);
// Average of a vector in the range [start, end)
double Average(const vector<double>& X, int start, int end);
// Compute probability density at position y
double Psi2(double y,double mu,double sigma);
// Compute local energy at position v
double LocalEnergy(double y, double mu, double sigma);

int main (int argc, char *argv[]) {

// Initialize the random number generator

	Random random;
	int seed[4];
	int p1, p2;

	// Read the prime addend (c) used to select the random stream
	ifstream inPrime("../Input/Primes");
	if (!inPrime.is_open()) {
		cerr<<"Input error: Unable to open Primes"<<endl;
		return 1;
	}
	inPrime>>p1>>p2;
	inPrime.close();

	// Read the initial seed and set the generator state
	ifstream inSeed("../Input/seed.in");
	if (!inSeed.is_open()) {
		cerr<<"Input error: Unable to open seed.in"<<endl;
		return 1;
	}
	string Property;
	while (inSeed>>Property) {
		if (Property=="RANDOMSEED") {
			inSeed>>seed[0]>>seed[1]>>seed[2]>>seed[3];
			random.SetRandom(seed,p1,p2);
		}
	}
	inSeed.close();	


	// Simulation parameters
	// Metropolis algorithm
	bool restart=false;						// If true, restart from last position
	double mu=-0.860813, sigma=0.582034;		// Wavefunction parameters: mean value, standard deviation of the two Gaussian distributions
	double x0=0.0;							// Starting position
	double delta=1.2;						// Metropolis step size
	int M=1000000;							// Number of total throws
	int N=1000;								// Number of blocks
	int L=M/N;								// Number of throws per block (block length)
	double y, x;							// Actual and proposed position
	int accepted=0;							// Number of accepted moves
	int attempted=0;						// Number of attempted moves

	// Simulated Annealing parameters
	bool annealing=false;					// If true, do Simulated Annealing
	double beta=1;							// Initial value of beta=1/T
	double dbeta=1;						// Beta increment
	int steps_beta=50;						// Number of optimization steps for each beta value
	double delta_mu=0.5, delta_sigma=0.5;	// Mu and sigma step sizes
	double err_mu=0.01, err_sigma=0.01;		// Required precision on mu and sigma

	// Arrays for data storage and blocking
	vector<double> HBuffer(L);				// Temporary buffer for energy values within one block
	vector<double> HBlock(N);				// Block averages: HBlock[i] = average of block i (i=0-N-1)
	vector<double> HBlock2(N);				// Squares of block averages: HBlock2[i] = HBlock[i]^2
	vector<double> HMeanValue(N);			// Progressive averages: HMeanValue[i] = average estimated over the first (i+1) blocks
	vector<double> HStDevOfTheMean(N);		// Progressive statistical uncertainty: HDevOfTheMean[i] = uncertainty of HMeanValue[i]

// Input
    ifstream inChoice("../Input/input.dat");
    if (!inChoice.is_open()) {
		cerr<<"Input error: Unable to open input.dat"<<endl;
		return 1;
	}
	string key;
	int value;
	while (inChoice>>key) {
		if (key=="restart") {
			inChoice>>value;
			restart=static_cast<bool>(value);
		} else if (key=="annealing") {
			inChoice>>value;
			annealing=static_cast<bool>(value);
		}
	}
	inChoice.close();

// 08.1. Variational Monte Carlo with uniform Metropolis sampling: the probability density |\Psi_{\mu,\sigma}(x)|^2 is sampled using a uniform proposal distribution  centered at the current position.
// For fixed variational parameters (\mu,\sigma), the expectation value of the Hamiltonian <H> is computed from the generated sequence and estimated via data blocking.

	if (annealing==false) {

		// Starting position
		if (restart==false) {					// Don't restart from final position of previous execution
			y=x0;
		} else {								// Restart from final position of previous execution. Read initial position from position.dat
			ifstream inPosition("../Output/position.dat");
			if (!inPosition.is_open()) {
				cerr<<"Input error: Unable to open position.dat"<<endl;
				return 1;
			}
			string Line, FinalLine;
			while (getline(inPosition,Line)) if (!Line.empty()) FinalLine=Line;
			inPosition.close();
			istringstream iss(FinalLine);
			iss>>y;
		}

		// Output files
		ofstream out_energy, out_position;
		out_energy.open("../Output/energy.dat");
		if(!out_energy.is_open()) {
			cerr<<"Error: unable to open energy.dat"<<endl;
			return 1;
		}
		out_position.open("../Output/position.dat");
		if(!out_position.is_open()) {
			cerr<<"Error: unable to open position.dat"<<endl;
			return 1;
		}

		// Hamiltonian expectation value
		// Loop over blocks
		for (int i=0; i<N; i++) {
			// Loop over samples within block
			for (int j=0; j<L; j++) {
				MetropolisMove(accepted,attempted,x,y,delta,mu,sigma,random);
				HBuffer[j]=LocalEnergy(y,mu,sigma);
			}
			out_position<<y<<endl;
			// Compute block average
			HBlock[i]=Average(HBuffer,0,L);
			HBlock2[i]=HBlock[i]*HBlock[i];
			// Compute progressive averages and uncertainties
			if (i==0) {
				HMeanValue[i]=HBlock[i];
				HStDevOfTheMean[i]=0.0;										// After the first block the statistical uncertainty is not defined
			} else {
				HMeanValue[i]=Average(HBlock,0,i+1);						// Progressive average over the first (i+1) blocks
				HStDevOfTheMean[i]=sqrt((Average(HBlock2,0,i+1)-HMeanValue[i]*HMeanValue[i])/i);	// Statistical uncertainty of the mean, over the first (i+1) blocks
			}
		}

		// Save results to file
		out_energy<<"# Energy calculation - Metropolis algorithm\n";
		out_energy<<"# Total number of throws: M = "<<M<<"\n";
		out_energy<<"# Number of blocks: N = "<<N<<"\n";
		out_energy<<"# Throws per block: L = "<<L<<"\n";
		out_energy<<"# Initial position: x0 = "<<x0<<"\n";
		out_energy<<"# Uniform proposal half-width: delta = "<<delta<<"\n";
		out_energy<<"# Acceptance rate: A = accepted/attempted = "<<static_cast<double>(accepted)/attempted<<"\n";
		out_energy<<"# Wavefunction parameters: (mu,sigma) = ("<<mu<<", "<<sigma<<")\n";
		out_energy<<"# Columns:\n";
		out_energy<<"# 1) block index\n";
		out_energy<<"# 2) block average\n";
		out_energy<<"# 3) progressive average\n";
		out_energy<<"# 4) progressive statistical uncertainty\n";
		for (int i=0; i<N; i++) out_energy<<i+1<<" "<<HBlock[i]<<" "<<HMeanValue[i]<<" "<<HStDevOfTheMean[i]<<"\n";
		out_energy.close();
		out_position.close();
	}


// 08.2. Variational parameter optimization via Simulated Annealing: the variational parameters (\mu, \sigma) are treated as dynamical variables and updated using a Metropolis algorithm in parameter space.
// The variational energy <H> plays the role of the cost function, and an annealing schedule on \beta is used to progressively approach the minimum-energy configuration.

	else if (annealing==true) {
		
		// Output file
		ofstream out_annealing("../Output/annealing.dat");
		if(!out_annealing.is_open()) {
			cerr<<"Error: unable to open annealing.dat"<<endl;
			return 1;
		}
		out_annealing<<"# Simulated Annealing\n";
		out_annealing<<"# Total number of throws: M = "<<M<<"\n";
		out_annealing<<"# Number of blocks: N = "<<N<<"\n";
		out_annealing<<"# Throws per block: L = "<<L<<"\n";
		out_annealing<<"# Initial position: x0 = "<<x0<<"\n";
		out_annealing<<"# Uniform proposal half-width for Metropolis: delta = "<<delta<<"\n";
		out_annealing<<"# Initial parameters: (mu,sigma) = ("<<mu<<", "<<sigma<<")\n";
		out_annealing<<"# Initial value of beta = "<<beta<<"\n";
		out_annealing<<"# Beta increment = "<<dbeta<<"\n";
		out_annealing<<"# Number of optimization steps for each beta value = "<<steps_beta<<"\n";
		out_annealing<<"# delta_mu = "<<delta_mu<<"\n";
		out_annealing<<"# delta_sigma = "<<delta_sigma<<"\n";
		out_annealing<<"# Required precision on mu: err_mu = "<<err_mu<<"\n";
		out_annealing<<"# Required precision on sigma: err_sigma = "<<err_sigma<<"\n";
		out_annealing<<"# Columns:\n";
		out_annealing<<"# 1) beta\n";
		out_annealing<<"# 2) mu\n";
		out_annealing<<"# 3) sigma\n";
		out_annealing<<"# 4) H_old\n";
		out_annealing<<"# 4) err_old\n";
		out_annealing<<"# 5) Acceptance rate: A = accepted/attempted\n";

		// Initial variational energy
		double H_old, H_new;
		y=x0;
		for (int i=0; i<N; i++) {
			for (int j=0; j<L; j++) {
				MetropolisMove(accepted,attempted,x,y,delta,mu,sigma,random);
				HBuffer[j]=LocalEnergy(y,mu,sigma);
			}
			HBlock[i]=Average(HBuffer,0,L);
		}
		H_old=Average(HBlock,0,N);
		
		// Simulated annealing loop
		// Progressive cooling: increase beta
		double err_old=0.0, err_new=0.0;
		while (2*delta_mu/beta>err_mu || 2*delta_sigma/beta>err_sigma) {
			for (int k=0; k<steps_beta; k++) {
				// Propose a move in parameter space (mu, sigma)
				double dmu=random.Rannyu(-delta_mu,delta_mu)/beta;
				double dsigma=random.Rannyu(-delta_sigma,delta_sigma)/beta;
				double mu_new=mu+dmu;
				double sigma_new=sigma+dsigma;
				if (sigma_new<=0.0) continue;				// Enforce sigma > 0
				double mu_backup=mu;
				double sigma_backup=sigma;
				mu=mu_new;
				sigma=sigma_new;
				// Evaluate the variational energy for the proposed parameters
				y=x0;
				accepted=0;
				attempted=0;
				for (int i=0; i<N; i++) {
					for (int j=0; j<L; j++) {
						MetropolisMove(accepted,attempted,x,y,delta,mu,sigma,random);
						HBuffer[j]=LocalEnergy(y,mu,sigma);
					}
					HBlock[i]=Average(HBuffer,0,L);
				}
				H_new=Average(HBlock,0,N);					// New variational energy
				double H2_new=0.0;
				for (int i=0; i<N; i++) H2_new+=HBlock[i]*HBlock[i];
				H2_new/=N;
				err_new=sqrt((H2_new-H_new*H_new)/(N-1));	// Statistical uncertainty
				double dH=H_old-H_new;						// Energy difference	
				// Metropolis acceptance rule in parameter space
				if (random.Rannyu()<Min(1.0,exp(beta*dH))) {
					H_old=H_new;
					err_old=err_new;
				}
				else {
					mu=mu_backup;							// Reject: restore previous parameters
					sigma=sigma_backup;
				}
				out_annealing<<beta<<" "<<mu<<" "<<sigma<<" "<<H_old<<" "<<err_old<<" "<<static_cast<double>(accepted)/attempted<<endl;
			}
			beta+=dbeta;									// Increase beta
		}

		ofstream out_af("../Output/annealing_final.dat");
		if(!out_af.is_open()) {
			cerr<<"Error: unable to open annealing_final.dat"<<endl;
			return 1;
		}
		out_af<<"Optimization results:\n";
		out_af<<"H = "<<H_old<<"\n";
		out_af<<"mu = "<<mu<<"\n";
		out_af<<"sigma = "<<sigma<<"\n";
		out_annealing.close();
		out_af.close();
	}

// Close the program

	random.SaveSeed();
	return 0; 

}

void MetropolisMove(int& accepted, int& attempted, double& x, double& y, double delta, double mu, double sigma, Random& random) {
	x=y+random.Rannyu(-delta,delta);						// New position proposal
	double alpha=Min(1.0,Psi2(x,mu,sigma)/Psi2(y,mu,sigma));
	if (random.Rannyu()<alpha) {			// Accept/reject
		y=x;
		accepted++;
	}
	attempted++;
}

double Min(double a, double b) {
	return (a<b)? a:b;
}

double Average(const vector<double>& X, int start, int end) {
	if (start<0||end<=start||end>static_cast<int>(X.size())) {
		cerr<<"Average error: invalid start/end indices"<<endl;
		return 0.0;
	}
	double sum=0.0;
	for (int i=start; i<end; i++) sum+=X[i];
	return sum/(end-start);
}

double Psi2(double y,double mu,double sigma) {
	double Psi=exp(-0.5*pow(y-mu,2)/pow(sigma,2))+exp(-0.5*pow(y+mu,2)/pow(sigma,2));
	return Psi*Psi;
}


double LocalEnergy(double y,double mu,double sigma) {
	double e1=exp(-0.5*pow(y-mu,2)/pow(sigma,2));
	double e2=exp(-0.5*pow(y+mu,2)/pow(sigma,2));
	double Psi=e1+e2;
	double K=-0.5*(e1*(pow(y-mu,2)/pow(sigma,4)-1.0/(sigma*sigma))+e2*(pow(y+mu,2)/pow(sigma,4)-1.0/(sigma*sigma)));
	double U=pow(y,4)-2.5*pow(y,2);
	return K/Psi + U;
}