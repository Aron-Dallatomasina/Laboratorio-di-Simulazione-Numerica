/********************************************************************************************************************************
 * Numerical Simulation Laboratory - Exercise 01.3
 * Monte Carlo simulation of Buffon's needle experiment for the estimation of π
 *
 * The program simulates the random throwing of a needle between parallel lines.
 * The probability of the needle intersecting one of the lines is estimated and used to compute π.
 * Statistical uncertainties are evaluated using the data blocking technique.
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

using namespace std;

// Average of a vector in the range [start, end)
double Average(const vector<double>& X, int start, int end);


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

// Simulation of the Buffon's experiment: generate needle throws and estimate <pi> using the blocking method

	// Simulation parameters
	int M=10000000;							// Total number of needle throws
	int N=1000;								// Number of blocks. Each block yields a measure of pi
	int N_throws=M/N;						// Number of throws per block (block length)
	double L=1;								// Needle length
	double d=1.2;							// Distance between parallel lines (d>L)
	// Arrays for data blocking
	vector<double> piBlock(N);				// pi estimate in each block: piBlock[i] = pi estimate in block i (i=0-N-1)
	vector<double> piBlock2(N);				// Squared block estimates: piBlock2[i] = piBlock[i]^2
	vector<double> piMeanValue(N);			// Progressive averages <pi>: piMeanValue[i] = <pi> estimated over the first (i+1) blocks
	vector<double> piStDevOfTheMean(N);		// Progressive statistical uncertainty of <pi>: piStDevOfTheMean[i] = uncertainty of piMeanValue[i]
    // Variables used during the simulation
	int N_hit;								// Number of hits in a block
	double h;								// y-coordinate of the needle center (uniformly distributed in [-d/2,d/2))
	double x, y;							// Random point used to define the needle direction (rejection sampling)
	double r2;								// x^2 + y^2

	// Simulation
	for (int i=0; i<N; i++) {
		N_hit=0;
		for (int j=0; j<N_throws; j++) {							// Block
			h=random.Rannyu(-0.5*d,0.5*d);							// Random position of the needle's center along y between two lines
			do {													// Generate a random direction without using pi: sample a point uniformly in the right half of the unit disk (x>=0) by rejection
				x=random.Rannyu(0.0,1.0);
				y=random.Rannyu(-1.0,1.0);
				r2=x*x+y*y;
			} while (r2>1.0||r2==0.0);
			if ( (h+fabs(y)/sqrt(r2)*L/2.0>d/2.0) || (h-fabs(y)/sqrt(r2)*L/2.0<-d/2.0) ) N_hit++;	// Needle intersects a line if one endpoint crosses y = +/-d/2
		}
		// Block estimate of pi from P=N_hit/N_throws and pi=2L/(dP)
		if (N_hit>0) piBlock[i]=2.0*L*N_throws/(d*N_hit);
		else piBlock[i]=0.0;
		// Data blocking: progressive mean and uncertainty
		piBlock2[i]=piBlock[i]*piBlock[i];
		if (i==0) {
			piMeanValue[i]=piBlock[i];
			piStDevOfTheMean[i]=0.0;								// After the first block the statistical uncertainty is not defined
		} else {
			piMeanValue[i]=Average(piBlock,0,i+1);					// Progressive average over the first (i+1) blocks
			piStDevOfTheMean[i]=sqrt((Average(piBlock2,0,i+1)-piMeanValue[i]*piMeanValue[i])/i);	// Statistical uncertainty of the mean, over the first (i+1) blocks
		}
	}

// Save results to file

	ofstream out("../Output/buffon_pi.dat");
	if (!out.is_open()) {
		cerr<<"Error: unable to open buffon_pi.dat"<<endl;
		return 1;
	}
	out<<"# Buffon's needle experiment - Monte Carlo estimation of pi\n";
	out<<"# Total number of throws: M = "<<M<<"\n";
	out<<"# Number of blocks: N = "<<N<<"\n";
	out<<"# Throws per block: N_throws = "<<N_throws<<"\n";
	out<<"# Needle length: L = "<<L<<"\n";
	out<<"# Distance between lines: d = "<<d<<"\n";
	out<<"# Columns:\n";
	out<<"# 1) block index\n";
	out<<"# 2) pi estimate in the block\n";
	out<<"# 3) progressive mean of pi\n";
	out<<"# 4) progressive statistical uncertainty\n\n";

	for (int i=0; i<N; i++) out<<i+1<<" "<<piBlock[i]<<" "<<piMeanValue[i]<<" "<<piStDevOfTheMean[i]<<"\n";
	out.close();

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
