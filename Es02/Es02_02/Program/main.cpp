/********************************************************************************************************************************
 * Numerical Simulation Laboratory - Exercise 02.2
 * Monte Carlo simulation of three-dimensional random walks
 *
 * The program simulates three-dimensional random walks always starting from the origin.
 * Two different models are considered:
 * 	1. Random walk on a cubic lattice with lattice constant a = 1
 * 	2. Random walk in the continuum with step length a = 1 and random directions uniformly distributed over the solid angle
 * The simulation is repeated many times and the root mean square distance from the origin sqrt(<|r_k|^2>) is evaluated as a function of the step k using the data blocking technique.
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
// Perform a single step of a random walk on a cubic lattice with lattice constant a=step. The walker moves by ±step along one of the three principal lattice directions (x,y,z)
void LatticeStep(Random &random, double &x, double &y, double &z, double step);

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

// 1. Simulation of a random walk on a cubic lattice with lattice constant a = 1:
// the walker starts from the origin and, at each step, moves by ±a along one of the three principal lattice directions (x, y, z). The simulation is repeated many times and the root mean square distance from the origin sqrt(<|r_k|^2>) is evaluated as a function of the step k.

	// Simulation parameters
	int M=1000000;							// Total number of throws (random walks)
	int N_steps=100;						// Number of steps for each random walk
	int N=1000;								// Number of blocks. Each block yields an estimate of the distance sqrt(<|r_k|^2>), for each step k
	int L=M/N;								// Number of throws per block (block length)
	double x, y, z;							// Cartesian coordinates of the walker
	double a=1.0;							// Step length (lattice constant)

	// Arrays for data storage and blocking
	vector<vector<double>> distanceBuffer(N_steps+1,vector<double>(L));		// Temporary buffer for values (squared distances) generated in one block - row index=step index, column index=throw (random walk) index. j-th column stores squared distance values for j-th random walk, for each step 
	vector<vector<double>> rmsdistanceBlock(N_steps+1,vector<double>(N));	// Square roots of block averages - row index=step index, column index=block index. rmsdistanceBlock[k][i] = square root of the average squared distance at k-th step, for block i (i=0-N-1) 
	vector<vector<double>> rmsdistanceBlock2(N_steps+1,vector<double>(N));	// rmsdistanceBlock2[k][i] = rmsdistanceBlock[k][i]^2
	vector<double> rmsdistanceMeanValue(N_steps+1);							// Averages of sqrt(<|r_k|^2>) over blocks, for each step k: rmsdistanceMeanValue[k] = root-mean-squared distance of the walker after k steps, estimated with the blocking method
	vector<double> rmsdistanceStDevOfTheMean(N_steps+1);					// Statistical uncertainties of < sqrt(<|r_k|^2>) > for each step k: rmsdistanceStDevOfTheMean[k] = uncertainty of rmsdistanceMeanValue[k]

	// Simulation
	for (int i=0; i<N; i++) {															// Loop over blocks

		for (int j=0; j<L; j++) {														// Loop over random walks in a block
			x=y=z=0.0;																	// Reset coordinates before starting the random walk
			for (int k=1; k<=N_steps; k++) {											// Loop over steps of a random walk
				LatticeStep(random,x,y,z,a);											// Take the k-th step
				distanceBuffer[k][j]=x*x+y*y+z*z;										// Compute squared distance after k-th step
			}
		}
		for (int k=1; k<=N_steps; k++) {												// Loop over steps
			rmsdistanceBlock[k][i]=sqrt(Average(distanceBuffer[k],0,L));				// Compute rms distance after k-th step for all the random walks in the i-th block
			rmsdistanceBlock2[k][i]=rmsdistanceBlock[k][i]*rmsdistanceBlock[k][i];		// Compute its squared value - used in computing statistical uncertainty
		}
	}

	rmsdistanceMeanValue[0]=0.0;														// After 0 steps, distance=0
	rmsdistanceStDevOfTheMean[0]=0.0;													// Same for uncertainty
	for (int k=1; k<=N_steps; k++) {
		rmsdistanceMeanValue[k]=Average(rmsdistanceBlock[k],0,N);						// Average rms distance after k-th step (averaging over blocks)
		rmsdistanceStDevOfTheMean[k]=sqrt((Average(rmsdistanceBlock2[k],0,N)-rmsdistanceMeanValue[k]*rmsdistanceMeanValue[k])/(N-1));	// Statistical uncertainty of the mean, after k-th step
	}

// Save results to file

	ofstream out_lattice("../Output/results_lattice.dat");
	if (!out_lattice.is_open()) {
		cerr<<"Error: unable to open results_lattice.dat"<<endl;
		return 1;
	}
	out_lattice<<"# Monte Carlo simulation of a 3D random walk on a cubic lattice\n";
	out_lattice<<"# Total number of throws (random walks): M = "<<M<<"\n";
	out_lattice<<"# Number of blocks: N = "<<N<<"\n";
	out_lattice<<"# Throws per block: L = "<<L<<"\n";
	out_lattice<<"# Number of steps: N_steps = "<<N_steps<<"\n";
	out_lattice<<"# Step length: a = "<<a<<"\n";
	out_lattice<<"# Columns:\n";
	out_lattice<<"# 1) step index\n";
	out_lattice<<"# 2) root mean square distance sqrt(<|r_k|^2>)\n";
	out_lattice<<"# 3) statistical uncertainty on sqrt(<|r_k|^2>)\n";

	for (int i=0; i<=N_steps; i++) out_lattice<<i<<" "<<rmsdistanceMeanValue[i]<<" "<<rmsdistanceStDevOfTheMean[i]<<"\n";
	out_lattice.close();

// 2. Simulation of a random walk in the continuum with step length a = 1:
// the walker starts from the origin and, at each step, moves by a along a random direction uniformly distributed over the solid angle. The simulation is repeated many times and the root mean square distance from the origin sqrt(<|r_k|^2>) is evaluated as a function of the step k.

	// Simulation parameters
	M=1000000;								// Total number of throws (random walks)
	N_steps=100;							// Number of steps for each random walk
	N=1000;									// Number of blocks. Each block yields an estimate of the distance sqrt(<|r_k|^2>), for each step k
	L=M/N;									// Number of throws per block (block length)
	x=y=z=0.0;								// Cartesian coordinates of the walker
	a=1.0;									// Step length

	// Arrays for data storage and blocking
	for (size_t i=0;i<distanceBuffer.size();i++) distanceBuffer[i].assign(L,0.0);			// Temporary buffer for values (squared distances) generated in one block
	for (size_t i=0;i<rmsdistanceBlock.size();i++) rmsdistanceBlock[i].assign(N,0.0);		// Square roots of block averages
	for (size_t i=0;i<rmsdistanceBlock2.size();i++) rmsdistanceBlock2[i].assign(N,0.0);	// rmsdistanceBlock2[k][i] = rmsdistanceBlock[k][i]^2
	rmsdistanceMeanValue.assign(N_steps+1,0.0);											// Averages of sqrt(<|r_k|^2>) over blocks, for each step k
	rmsdistanceStDevOfTheMean.assign(N_steps+1,0.0);									// Statistical uncertainties of < sqrt(<|r_k|^2>) > for each step k

	// Simulation
	for (int i=0; i<N; i++) {															// Loop over blocks

		for (int j=0; j<L; j++) {														// Loop over random walks in a block
			x=y=z=0.0;																	// Reset coordinates before starting the random walk
			for (int k=1; k<=N_steps; k++) {											// Loop over steps of a random walk
				double theta=random.Theta();											// Take the k-th step
				double phi=random.Rannyu(0,2.0*pi);
				x+=a*sin(theta)*cos(phi);
				y+=a*sin(theta)*sin(phi);
				z+=a*cos(theta);
				distanceBuffer[k][j]=x*x+y*y+z*z;										// Compute squared distance after k-th step
			}
		}
		for (int k=1; k<=N_steps; k++) {												// Loop over steps
			rmsdistanceBlock[k][i]=sqrt(Average(distanceBuffer[k],0,L));				// Compute rms distance after k-th step for all the random walks in the i-th block
			rmsdistanceBlock2[k][i]=rmsdistanceBlock[k][i]*rmsdistanceBlock[k][i];		// Compute its squared value - used in computing statistical uncertainty
		}
	}

	rmsdistanceMeanValue[0]=0.0;														// After 0 steps, distance=0
	rmsdistanceStDevOfTheMean[0]=0.0;													// Same for uncertainty
	for (int k=1; k<=N_steps; k++) {
		rmsdistanceMeanValue[k]=Average(rmsdistanceBlock[k],0,N);						// Average rms distance after k-th step (averaging over blocks)
		rmsdistanceStDevOfTheMean[k]=sqrt((Average(rmsdistanceBlock2[k],0,N)-rmsdistanceMeanValue[k]*rmsdistanceMeanValue[k])/(N-1));	// Statistical uncertainty of the mean, after k-th step
	}

// Save results to file

	ofstream out_continuum("../Output/results_continuum.dat");
	if (!out_continuum.is_open()) {
		cerr<<"Error: unable to open results_continuum.dat"<<endl;
		return 1;
	}
	out_continuum<<"# Monte Carlo simulation of a 3D random walk in the continuum\n";
	out_continuum<<"# Total number of throws (random walks): M = "<<M<<"\n";
	out_continuum<<"# Number of blocks: N = "<<N<<"\n";
	out_continuum<<"# Throws per block: L = "<<L<<"\n";
	out_continuum<<"# Number of steps: N_steps = "<<N_steps<<"\n";
	out_continuum<<"# Step length: a = "<<a<<"\n";
	out_continuum<<"# Columns:\n";
	out_continuum<<"# 1) step index\n";
	out_continuum<<"# 2) root mean square distance sqrt(<|r_k|^2>)\n";
	out_continuum<<"# 3) statistical uncertainty on sqrt(<|r_k|^2>)\n";

	for (int i=0; i<=N_steps; i++) out_continuum<<i<<" "<<rmsdistanceMeanValue[i]<<" "<<rmsdistanceStDevOfTheMean[i]<<"\n";
	out_continuum.close();

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

void LatticeStep(Random &random, double &x, double &y, double &z, double step) {
	int s=random.UniformDice();
	if (s==1)		x+=step;
	else if (s==2)	x-=step;
	else if (s==3)	y+=step;
	else if (s==4)	y-=step;
	else if (s==5)	z+=step;
	else if (s==6)	z-=step;
}