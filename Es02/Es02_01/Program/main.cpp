/********************************************************************************************************************************
 * Numerical Simulation Laboratory - Exercise 02.1
 * Monte Carlo estimation of a one-dimensional integral with uniform sampling and importance sampling
 *
 * The program estimates the integral I = ∫_0^1 (π/2) cos(πx/2) dx by Monte Carlo integration using
 * 	1. Uniform sampling in [0,1)
 * 	2. Importance sampling using a non-uniform probability density in [0,1)
 * For both methods, progressive averages and statistical uncertainties are evaluated with the data blocking technique.
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
// Integrand function f(x)=(π/2)cos(πx/2) for the Monte Carlo integration
double f(double x);
// Auxiliary function g(x)=f(x)/ρ(x)=cos(πx/2)/(2/π+1/2-x) used for importance sampling
double g(double x);

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

// 1. Monte Carlo estimation of the integral (I) using uniform sampling in [0,1): the integrand is evaluated on uniformly distributed random points and the integral is computed as the average of the sampled values

	// Simulation parameters
	int M=10000000;							// Total number of throws (integral evaluations)
	int N=1000;								// Number of blocks. Each block yields a measure of I
	int L=M/N;								// Number of throws per block (block length)
	// Arrays for data blocking
	vector<double> IBuffer(L);				// Temporary buffer for values generated in one block
	vector<double> IBlock(N);				// Block averages: IBlock[i] = average of block i (i=0-N-1)
	vector<double> IBlock2(N);				// Squares of block averages: IBlock2[i] = IBlock[i]^2
	vector<double> IMeanValue(N);			// Progressive averages <I>: IMeanValue[i] = <I> estimated over the first (i+1) blocks
	vector<double> IStDevOfTheMean(N);		// Progressive statistical uncertainty of <I>: IStDevOfTheMean[i] = uncertainty of IMeanValue[i]

	// Simulation
	for (int i=0; i<N; i++) {
		for (int j=0; j<L; j++) IBuffer[j]=f(random.Rannyu());		// Generate L random numbers for the i-th block and evaluate the integrand in all of them
		IBlock[i]=Average(IBuffer,0,L);								// Compute the average of the i-th block
		IBlock2[i]=IBlock[i]*IBlock[i];
		if (i==0) {
			IMeanValue[i]=IBlock[i];
			IStDevOfTheMean[i]=0.0;									// After the first block the statistical uncertainty is not defined
		} else {
			IMeanValue[i]=Average(IBlock,0,i+1);					// Progressive average over the first (i+1) blocks
			IStDevOfTheMean[i]=sqrt((Average(IBlock2,0,i+1)-IMeanValue[i]*IMeanValue[i])/i);	// Statistical uncertainty of the mean, over the first (i+1) blocks
		}
	}

// Save results to file

	ofstream out_uniform("../Output/I_uniform_sampling.dat");
	if (!out_uniform.is_open()) {
		cerr<<"Error: unable to open I_uniform_sampling.dat"<<endl;
		return 1;
	}
	out_uniform<<"# Monte Carlo estimation of the integral using uniform sampling in [0,1)\n";
	out_uniform<<"# Total number of throws: M = "<<M<<"\n";
	out_uniform<<"# Number of blocks: N = "<<N<<"\n";
	out_uniform<<"# Throws per block: L = "<<L<<"\n";
	out_uniform<<"# Columns:\n";
	out_uniform<<"# 1) block index\n";
	out_uniform<<"# 2) I estimate in the block\n";
	out_uniform<<"# 3) progressive mean of I\n";
	out_uniform<<"# 4) progressive statistical uncertainty\n\n";

	for (int i=0; i<N; i++) out_uniform<<i+1<<" "<<IBlock[i]<<" "<<IMeanValue[i]<<" "<<IStDevOfTheMean[i]<<"\n";
	out_uniform.close();

// 2. Monte Carlo estimation of the integral (I) using importance sampling in [0,1): the integrand is written as f(x)=g(x)ρ(x), and the integral is computed as the average of g(x) over random points distributed according to ρ(x). The chosen pdf is a linear approximation of the integrand obtained from its first-order Taylor expansion around x=1/2: ρ(x)=(1+π/4)-(π/2)*x

	random.SetRandom(seed,p1,p2);
	// Simulation parameters
	M=10000000;								// Total number of throws (integral evaluations)
	N=1000;									// Number of blocks. Each block yields a measure of I
	L=M/N;									// Number of throws per block (block length)
	// Arrays for data blocking
	IBuffer.assign(L,0.0);					// Temporary buffer for values generated in one block
	IBlock.assign(N,0.0);					// Block averages: IBlock[i] = average of block i (i=0-N-1)
	IBlock2.assign(N,0.0);					// Squares of block averages: IBlock2[i] = IBlock[i]^2
	IMeanValue.assign(N,0.0);				// Progressive averages <I>: IMeanValue[i] = <I> estimated over the first (i+1) blocks
	IStDevOfTheMean.assign(N,0.0);			// Progressive statistical uncertainty of <I>: IStDevOfTheMean[i] = uncertainty of IMeanValue[i]

	// Simulation
	double m=-pi/2.0;						// Angular coefficient of the chosen pdf ρ(x): m=-π/2
	for (int i=0; i<N; i++) {
		for (int j=0; j<L; j++) IBuffer[j]=g(random.LinearIn01(m));	// Generate L random numbers for the i-th block and evaluate the auxiliary function g(x) in all of them
		IBlock[i]=Average(IBuffer,0,L);								// Compute the average of the i-th block
		IBlock2[i]=IBlock[i]*IBlock[i];
		if (i==0) {
			IMeanValue[i]=IBlock[i];
			IStDevOfTheMean[i]=0.0;									// After the first block the statistical uncertainty is not defined
		} else {
			IMeanValue[i]=Average(IBlock,0,i+1);					// Progressive average over the first (i+1) blocks
			IStDevOfTheMean[i]=sqrt((Average(IBlock2,0,i+1)-IMeanValue[i]*IMeanValue[i])/i);	// Statistical uncertainty of the mean, over the first (i+1) blocks
		}
	}

// Save results to file

	ofstream out_importance("../Output/I_importance_sampling.dat");
	if (!out_importance.is_open()) {
		cerr<<"Error: unable to open I_importance_sampling.dat"<<endl;
		return 1;
	}
	out_importance<<"# Monte Carlo estimation of the integral using importance sampling with pdf (1+π/4)-(π/2)*x\n";
	out_importance<<"# Total number of throws: M = "<<M<<"\n";
	out_importance<<"# Number of blocks: N = "<<N<<"\n";
	out_importance<<"# Throws per block: L = "<<L<<"\n";
	out_importance<<"# Columns:\n";
	out_importance<<"# 1) block index\n";
	out_importance<<"# 2) I estimate in the block\n";
	out_importance<<"# 3) progressive mean of I\n";
	out_importance<<"# 4) progressive statistical uncertainty\n\n";

	for (int i=0; i<N; i++) out_importance<<i+1<<" "<<IBlock[i]<<" "<<IMeanValue[i]<<" "<<IStDevOfTheMean[i]<<"\n";
	out_importance.close();

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

double f(double x) {
	return pi/2.0*cos(pi*x/2.0);
}

double g(double x) {
	return cos(pi*x/2.0) / (2.0/pi+0.5-x);
}