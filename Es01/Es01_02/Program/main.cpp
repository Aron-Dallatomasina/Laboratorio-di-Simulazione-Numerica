/********************************************************************************************************************************
 * Numerical Simulation Laboratory - Exercise 01.2
 * Extension of the RANNYU pseudo-random number generator and verification of the Central Limit Theorem
 *
 * The program implements exponential and Cauchy–Lorentz distributions using the inversion of the cumulative distribution function.
 * It then generates data to study the behavior of sample averages and to verify the Central Limit Theorem for different probability distributions.
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

// Compute the averages S_N of N = 1, 2, 10, 100 random numbers extracted from different distributions

	vector<int> N{1,2,10,100};				// Number of random variables used to compute each average
	int M=10000;							// Number of averages to be generated for each value of N
	double sum;								// Temporary variable used to accumulate the sum

	// Uniform distribution in [0,1), generated with Random::Rannyu()
	vector<vector<double>> MeanRannyu(int(N.size()),vector<double>(M));				// MeanRannyu[n][i] = i-th realization of S_N for N=N[n]
	for (int n=0; n<int(N.size()); n++) {
		for (int i=0; i<M; i++) {
			sum=0.0;
			for (int j=0; j<N[n]; j++) sum+=random.Rannyu();
			MeanRannyu[n][i]=sum/N[n];
		}
	}

	// Exponential distribution with lambda=1, generated with Random::Exponential()
	vector<vector<double>> MeanExponential(int(N.size()),vector<double>(M));		// MeanExponential[n][i] = i-th realization of S_N for N=N[n]
	for (int n=0; n<int(N.size()); n++) {
		for (int i=0; i<M; i++) {
			sum=0.0;
			for (int j=0; j<N[n]; j++) sum+=random.Exponential(1.0);
			MeanExponential[n][i]=sum/N[n];
		}
	}

	// Cauchy–Lorentz distribution with mu=0 and Gamma=1, generated with Random::CauchyLorentz()
	vector<vector<double>> MeanCauchyLorentz(int(N.size()),vector<double>(M));		// MeanCauchyLorentz[n][i] = i-th realization of S_N for N=N[n]
	for (int n=0; n<int(N.size()); n++) {
		for (int i=0; i<M; i++) {
			sum=0.0;
			for (int j=0; j<N[n]; j++) sum+=random.CauchyLorentz(0.0,1.0);
			MeanCauchyLorentz[n][i]=sum/N[n];
		}
	}

// Save results to file

	// Uniform distribution
	ofstream out_rannyu("../Output/CLT_rannyu.dat");
	if (!out_rannyu) {
		cerr<<"Output error: Unable to open CLT_rannyu.dat"<<endl;
		return 1;
	}
	out_rannyu<<"# Central Limit Theorem - Uniform distribution in [0,1)\n";
	out_rannyu<<"# Columns: S_N for N = 1, 2, 10, 100\n";
	out_rannyu<<"# Number of realizations per N: M = "<<M<<"\n";
	for (int i=0; i<M; i++) out_rannyu<<MeanRannyu[0][i]<<" "<<MeanRannyu[1][i]<<" "<<MeanRannyu[2][i]<<" "<<MeanRannyu[3][i]<<"\n";
	out_rannyu.close();

	// Exponential distribution
	ofstream out_exp("../Output/CLT_exponential.dat");
	if (!out_exp) {
		cerr<<"Output error: Unable to open CLT_exponential.dat"<<endl;
		return 1;
	}
	out_exp<<"# Central Limit Theorem - Exponential distribution (lambda=1)\n";
	out_exp<<"# Columns: S_N for N = 1, 2, 10, 100\n";
	out_exp<<"# Number of realizations per N: M = "<<M<<"\n";
	for (int i=0; i<M; i++) out_exp<<MeanExponential[0][i]<<" "<<MeanExponential[1][i]<<" "<<MeanExponential[2][i]<<" "<<MeanExponential[3][i]<<"\n";
	out_exp.close();

	// Cauchy–Lorentz distribution
	ofstream out_cauchy("../Output/CLT_cauchy.dat");
	if (!out_cauchy) {
		cerr<<"Output error: Unable to open CLT_cauchy.dat"<<endl;
		return 1;
	}
	out_cauchy<<"# Central Limit Theorem - Cauchy-Lorentz distribution (mu=0, Gamma=1)\n";
	out_cauchy<<"# Columns: S_N for N = 1, 2, 10, 100\n";
	out_cauchy<<"# Number of realizations per N: M = "<<M<<"\n";
	for (int i=0; i<M; i++) out_cauchy<<MeanCauchyLorentz[0][i]<<" "<<MeanCauchyLorentz[1][i]<<" "<<MeanCauchyLorentz[2][i]<<" "<<MeanCauchyLorentz[3][i]<<"\n";
	out_cauchy.close();

// Close the program

	random.SaveSeed();
	return 0; 

}

