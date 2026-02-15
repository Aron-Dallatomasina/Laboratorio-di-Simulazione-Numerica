/********************************************************************************************************************************
 * Numerical Simulation Laboratory - Exercise 05.1
 * Metropolis sampling of hydrogen atom wavefunctions
 *
 * The program samples the probability densities |\Psi(x,y,z)|^2 of the hydrogen atom (1s and 2p states) using the Metropolis algorithm in three dimensions.
 *
 * Uniform and Gaussian proposal distributions are implemented. The expectation value <r> is estimated via data blocking.
 ********************************************************************************************************************************/



 // C++ Standard Library
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <iomanip>
#include <cmath>
#include <string>
// C++ Standard Template Library
#include <vector>
#include <algorithm>
#include <numeric>
// RNG
#include "random.h"

constexpr double pi=3.14159265358979323846;
using namespace std;

// Average of a vector in the range [start, end)
double Average(const vector<double>& X, int start, int end);
// Return the minimum between two real numbers a and b
double Min(double a, double b);
// Return the ratio rho(trial)/rho(current). 1s and 2p cases
double Ratio1s(double x_trial, double y_trial, double z_trial, double x_current, double y_current, double z_current);
double Ratio2p(double x_trial, double y_trial, double z_trial, double x_current, double y_current, double z_current);



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

// 1. Metropolis sampling with uniform transition probability:
// the probability densities |\Psi_{1,0,0}(x,y,z)|^2 and |\Psi_{2,1,0}(x,y,z)|^2 are sampled using a uniform proposal distribution in a cube centered at the current point. The expectation value <r> is computed from the generated sequence and estimated via data blocking.

	// Simulation parameters
	int M=10000000;							// Total number of sampled points
	int N=10000;							// Number of blocks. Each block yields an estimate of <r>
	int L=M/N;								// Number of samples per block (block length)
	double a_1s=1.20;						// Half side length of the uniform proposal cube (reduced units) for 1s state
	double a_2p=3.00;						// Half side length of the uniform proposal cube (reduced units) for 2p state
	double x_01s=0.0;						// Initial x-coordinate for 1s state (reduced units)
	double y_01s=0.0;						// y
	double z_01s=0.0;						// z
	double x_02p=0.0;						// Initial x-coordinate for 2p state (reduced units)
	double y_02p=0.0;						// y
	double z_02p=1.0;						// z

	// Arrays for data storage and blocking
	double x_old=0.0;						// Current x-coordinate
	double y_old=0.0;						// y
	double z_old=0.0;						// z
	double x_new=0.0;						// Trial x-coordinate
	double y_new=0.0;						// y
	double z_new=0.0;						// z
	int accepted_1s=0;						// Number of accepted moves (1s)
	int attempted_1s=0;						// Number of attempted moves (1s)
	int accepted_2p=0;						// Number of accepted moves (1s)
	int attempted_2p=0;						// Number of attempted moves (1s)
	// 1s state
	vector<double> r1sBuffer(L);			// Temporary buffer for radius values within one block (1s)
	vector<double> r1sBlock(N);				// Block averages <r> (1s): r1sBlock[i] = average of block i (i=0-N-1)
	vector<double> r1sBlock2(N);			// Squares of block averages (1s): r1sBlock2[i] = r1sBlock[i]^2
	vector<double> r1sMeanValue(N);			// Progressive averages of <r>, avg(<r>) (1s): r1sMeanValue[i] = avg(<r>) estimated over the first (i+1) blocks
	vector<double> r1sStDevOfTheMean(N);	// Progressive statistical uncertainty of avg(<r>) (1s): r1sDevOfTheMean[i] = uncertainty of r1sMeanValue[i]
	// 2p state
	vector<double> r2pBuffer(L);			// Temporary buffer for radius values within one block (2p)
	vector<double> r2pBlock(N);				// Block averages <r> (2p): r2pBlock[i] = average of block i (i=0-N-1)
	vector<double> r2pBlock2(N);			// Squares of block averages (2p): r2pBlock2[i] = r2pBlock[i]^2
	vector<double> r2pMeanValue(N);			// Progressive averages of <r>, avg(<r>) (2p): r2pMeanValue[i] = avg(<r>) estimated over the first (i+1) blocks
	vector<double> r2pStDevOfTheMean(N);	// Progressive statistical uncertainty of avg(<r>) (2p): r2pDevOfTheMean[i] = uncertainty of r2pMeanValue[i]

	// Sampling and computation of <r> for the 1s state
	// Initialize Markov chain at chosen starting point
	x_old=x_01s;
	y_old=y_01s;
	z_old=z_01s;
	// Loop over blocks
	for (int i=0; i<N; i++) {
		// Loop over samples within block
		for (int j=0; j<L; j++) {
			// Generate trial move from uniform proposal in cube of side 2*a_1s
			x_new=x_old+random.Rannyu(-a_1s,a_1s);
			y_new=y_old+random.Rannyu(-a_1s,a_1s);
			z_new=z_old+random.Rannyu(-a_1s,a_1s);
			attempted_1s++;
			// Metropolis acceptance step
			if (random.Rannyu()<=Min(1.0,Ratio1s(x_new,y_new,z_new,x_old,y_old,z_old))) {
				x_old=x_new;
				y_old=y_new;
				z_old=z_new;
				accepted_1s++;
			}
			// Store radius of current state
			r1sBuffer[j]=sqrt(x_old*x_old+y_old*y_old+z_old*z_old);
		}
		// Compute block average
		r1sBlock[i]=Average(r1sBuffer,0,L);
		r1sBlock2[i]=r1sBlock[i]*r1sBlock[i];
		// Compute progressive averages and uncertainties
		if (i==0) {
			r1sMeanValue[i]=r1sBlock[i];
			r1sStDevOfTheMean[i]=0.0;										// After the first block the statistical uncertainty is not defined
		} else {
			r1sMeanValue[i]=Average(r1sBlock,0,i+1);						// Progressive average over the first (i+1) blocks
			r1sStDevOfTheMean[i]=sqrt((Average(r1sBlock2,0,i+1)-r1sMeanValue[i]*r1sMeanValue[i])/i);	// Statistical uncertainty of the mean, over the first (i+1) blocks
		}

	}

	// Sampling and computation of <r> for the 2p state
	// Initialize Markov chain at chosen starting point
	x_old=x_02p;
	y_old=y_02p;
	z_old=z_02p;
	// Loop over blocks
	for (int i=0; i<N; i++) {
		// Loop over samples within block
		for (int j=0; j<L; j++) {
			// Generate trial move from uniform proposal in cube of side 2*a_2p
			x_new=x_old+random.Rannyu(-a_2p,a_2p);
			y_new=y_old+random.Rannyu(-a_2p,a_2p);
			z_new=z_old+random.Rannyu(-a_2p,a_2p);
			attempted_2p++;
			// Metropolis acceptance step
			if (random.Rannyu()<=Min(1.0,Ratio2p(x_new,y_new,z_new,x_old,y_old,z_old))) {
				x_old=x_new;
				y_old=y_new;
				z_old=z_new;
				accepted_2p++;
			}
			// Store radius of current state
			r2pBuffer[j]=sqrt(x_old*x_old+y_old*y_old+z_old*z_old);
		}
		// Compute block average
		r2pBlock[i]=Average(r2pBuffer,0,L);
		r2pBlock2[i]=r2pBlock[i]*r2pBlock[i];
		// Compute progressive averages and uncertainties
		if (i==0) {
			r2pMeanValue[i]=r2pBlock[i];
			r2pStDevOfTheMean[i]=0.0;										// After the first block the statistical uncertainty is not defined
		} else {
			r2pMeanValue[i]=Average(r2pBlock,0,i+1);						// Progressive average over the first (i+1) blocks
			r2pStDevOfTheMean[i]=sqrt((Average(r2pBlock2,0,i+1)-r2pMeanValue[i]*r2pMeanValue[i])/i);	// Statistical uncertainty of the mean, over the first (i+1) blocks
		}

	}

// Save results to file

	ofstream out_unif("../Output/results_unif.dat");
	if (!out_unif.is_open()) {
		cerr<<"Error: unable to open results_unif"<<endl;
		return 1;
	}
	out_unif<<"# Mean radius calculation for 1s and 2p (m=0) Hydrogen state - Metropolis algorithm with uniform transition probability\n";
	out_unif<<"# Total number of throws (sampled points): M = "<<M<<"\n";
	out_unif<<"# Number of blocks: N = "<<N<<"\n";
	out_unif<<"# Throws per block: L = "<<L<<"\n";
	out_unif<<"# 1s state ----------\n";
	out_unif<<"# Initial position: (x,y,z) = ("<<x_01s<<", "<<y_01s<<", "<<z_01s<<")\n";
	out_unif<<"# Uniform proposal half-width (cube side = 2a): a = "<<a_1s<<"\n";
	out_unif<<"# Acceptance rate: A = accepted/attempted = "<<static_cast<double>(accepted_1s)/attempted_1s<<"\n";
	out_unif<<"# 2p state ----------\n";
	out_unif<<"# Initial position: (x,y,z) = ("<<x_02p<<", "<<y_02p<<", "<<z_02p<<")\n";
	out_unif<<"# Uniform proposal half-width (cube side = 2a): a = "<<a_2p<<"\n";
	out_unif<<"# Acceptance rate: A = accepted/attempted = "<<static_cast<double>(accepted_2p)/attempted_2p<<"\n";
	out_unif<<"# Columns:\n";
	out_unif<<"# 1) block index\n";
	out_unif<<"# 2) block average for 1s state\n";
	out_unif<<"# 3) progressive average <r> for 1s state\n";
	out_unif<<"# 4) progressive statistical uncertainty for 1s state\n";
	out_unif<<"# 5) block average for 2p state\n";
	out_unif<<"# 6) progressive average <r> for 2p state\n";
	out_unif<<"# 7) progressive statistical uncertainty for 2p state\n";

	for (int i=0; i<N; i++) out_unif<<i+1<<" "<<r1sBlock[i]<<" "<<r1sMeanValue[i]<<" "<<r1sStDevOfTheMean[i]<<" "<<r2pBlock[i]<<" "<<r2pMeanValue[i]<<" "<<r2pStDevOfTheMean[i]<<"\n";
	out_unif.close();

// 2. Metropolis sampling with Gaussian transition probability:
// the probability densities |\Psi_{1,0,0}(x,y,z)|^2 and |\Psi_{2,1,0}(x,y,z)|^2 are sampled using a Gaussian proposal distribution centered at the current point. The expectation value <r> is computed from the generated sequence and estimated via data blocking.

// The following implementation is identical to the uniform case, except for
//	- the generation of the trial move (Gaussian instead of uniform)
//	- the proposal parameter (sigma instead of a).

	// Simulation parameters
	M=10000000;
	N=10000;
	L=M/N;
	double sigma_1s=0.75;						// Standard deviation of the Gaussian proposal distribution (reduced units) for 1s state
	double sigma_2p=1.90;						// Standard deviation of the Gaussian proposal distribution (reduced units) for 2p state
	x_01s=0.0;
	y_01s=0.0;
	z_01s=0.0;
	x_02p=0.0;
	y_02p=0.0;
	z_02p=1.0;

	// Arrays for data storage and blocking
	x_old=0.0;
	y_old=0.0;
	z_old=0.0;
	x_new=0.0;
	y_new=0.0;
	z_new=0.0;
	accepted_1s=attempted_1s=0;
	accepted_2p=attempted_2p=0;
	// 1s state
	r1sBuffer.assign(L,0.0);
	r1sBlock.assign(N,0.0);
	r1sBlock2.assign(N,0.0);
	r1sMeanValue.assign(N,0.0);
	r1sStDevOfTheMean.assign(N,0.0);
	// 2p state
	r2pBuffer.assign(L,0.0);
	r2pBlock.assign(N,0.0);
	r2pBlock2.assign(N,0.0);
	r2pMeanValue.assign(N,0.0);
	r2pStDevOfTheMean.assign(N,0.0);

	// Sampling and computation of <r> for the 1s state
	// Initialize Markov chain at chosen starting point
	x_old=x_01s;
	y_old=y_01s;
	z_old=z_01s;
	// Loop over blocks
	for (int i=0; i<N; i++) {
		// Loop over samples within block
		for (int j=0; j<L; j++) {
			// Generate trial move from Gaussian proposal centered at the current point with standard deviation sigma_1s
			x_new=x_old+random.Gauss(0.0,sigma_1s);
			y_new=y_old+random.Gauss(0.0,sigma_1s);
			z_new=z_old+random.Gauss(0.0,sigma_1s);
			attempted_1s++;
			// Metropolis acceptance step
			if (random.Rannyu()<=Min(1.0,Ratio1s(x_new,y_new,z_new,x_old,y_old,z_old))) {
				x_old=x_new;
				y_old=y_new;
				z_old=z_new;
				accepted_1s++;
			}
			// Store radius of current state
			r1sBuffer[j]=sqrt(x_old*x_old+y_old*y_old+z_old*z_old);
		}
		// Compute block average
		r1sBlock[i]=Average(r1sBuffer,0,L);
		r1sBlock2[i]=r1sBlock[i]*r1sBlock[i];
		// Compute progressive averages and uncertainties
		if (i==0) {
			r1sMeanValue[i]=r1sBlock[i];
			r1sStDevOfTheMean[i]=0.0;
		} else {
			r1sMeanValue[i]=Average(r1sBlock,0,i+1);
			r1sStDevOfTheMean[i]=sqrt((Average(r1sBlock2,0,i+1)-r1sMeanValue[i]*r1sMeanValue[i])/i);
		}

	}

	// Sampling and computation of <r> for the 2p state
	// Initialize Markov chain at chosen starting point
	x_old=x_02p;
	y_old=y_02p;
	z_old=z_02p;
	// Loop over blocks
	for (int i=0; i<N; i++) {
		// Loop over samples within block
		for (int j=0; j<L; j++) {
			// Generate trial move from Gaussian proposal centered at the current point with standard deviation sigma_2p
			x_new=x_old+random.Gauss(0.0,sigma_2p);
			y_new=y_old+random.Gauss(0.0,sigma_2p);
			z_new=z_old+random.Gauss(0.0,sigma_2p);
			attempted_2p++;
			// Metropolis acceptance step
			if (random.Rannyu()<=Min(1.0,Ratio2p(x_new,y_new,z_new,x_old,y_old,z_old))) {
				x_old=x_new;
				y_old=y_new;
				z_old=z_new;
				accepted_2p++;
			}
			// Store radius of current state
			r2pBuffer[j]=sqrt(x_old*x_old+y_old*y_old+z_old*z_old);
		}
		// Compute block average
		r2pBlock[i]=Average(r2pBuffer,0,L);
		r2pBlock2[i]=r2pBlock[i]*r2pBlock[i];
		// Compute progressive averages and uncertainties
		if (i==0) {
			r2pMeanValue[i]=r2pBlock[i];
			r2pStDevOfTheMean[i]=0.0;
		} else {
			r2pMeanValue[i]=Average(r2pBlock,0,i+1);
			r2pStDevOfTheMean[i]=sqrt((Average(r2pBlock2,0,i+1)-r2pMeanValue[i]*r2pMeanValue[i])/i);
		}

	}

// Save results to file

	ofstream out_gaus("../Output/results_gaus.dat");
	if (!out_gaus.is_open()) {
		cerr<<"Error: unable to open results_gaus"<<endl;
		return 1;
	}
	out_gaus<<"# Mean radius calculation for 1s and 2p (m=0) Hydrogen state - Metropolis algorithm with Gaussian transition probability\n";
	out_gaus<<"# Total number of throws (sampled points): M = "<<M<<"\n";
	out_gaus<<"# Number of blocks: N = "<<N<<"\n";
	out_gaus<<"# Throws per block: L = "<<L<<"\n";
	out_gaus<<"# 1s state ----------\n";
	out_gaus<<"# Initial position: (x,y,z) = ("<<x_01s<<", "<<y_01s<<", "<<z_01s<<")\n";
	out_gaus<<"# Gaussian proposal standard deviation: sigma = "<<sigma_1s<<"\n";
	out_gaus<<"# Acceptance rate: A = accepted/attempted = "<<static_cast<double>(accepted_1s)/attempted_1s<<"\n";
	out_gaus<<"# 2p state ----------\n";
	out_gaus<<"# Initial position: (x,y,z) = ("<<x_02p<<", "<<y_02p<<", "<<z_02p<<")\n";
	out_gaus<<"# Gaussian proposal standard deviation: sigma = "<<sigma_2p<<"\n";
	out_gaus<<"# Acceptance rate: A = accepted/attempted = "<<static_cast<double>(accepted_2p)/attempted_2p<<"\n";
	out_gaus<<"# Columns:\n";
	out_gaus<<"# 1) block index\n";
	out_gaus<<"# 2) block average for 1s state\n";
	out_gaus<<"# 3) progressive average <r> for 1s state\n";
	out_gaus<<"# 4) progressive statistical uncertainty for 1s state\n";
	out_gaus<<"# 5) block average for 2p state\n";
	out_gaus<<"# 6) progressive average <r> for 2p state\n";
	out_gaus<<"# 7) progressive statistical uncertainty for 2p state\n";

	for (int i=0; i<N; i++) out_gaus<<i+1<<" "<<r1sBlock[i]<<" "<<r1sMeanValue[i]<<" "<<r1sStDevOfTheMean[i]<<" "<<r2pBlock[i]<<" "<<r2pMeanValue[i]<<" "<<r2pStDevOfTheMean[i]<<"\n";
	out_gaus.close();


// Close the program

	random.SaveSeed();
	return 0; 

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

double Min(double a, double b) {
	return (a<b)? a:b;
}

double Ratio1s(double x_trial, double y_trial, double z_trial, double x_current, double y_current, double z_current) {
	return exp(2.0*(sqrt(x_current*x_current+y_current*y_current+z_current*z_current)-sqrt(x_trial*x_trial+y_trial*y_trial+z_trial*z_trial)));
}
double Ratio2p(double x_trial, double y_trial, double z_trial, double x_current, double y_current, double z_current) {
	if (fabs(z_current) < 1e-12) return 1.0;			// current point has rho ≈ 0 → accept any move
	return (z_trial/z_current)*(z_trial/z_current)*exp((sqrt(x_current*x_current+y_current*y_current+z_current*z_current)-sqrt(x_trial*x_trial+y_trial*y_trial+z_trial*z_trial)));
}
