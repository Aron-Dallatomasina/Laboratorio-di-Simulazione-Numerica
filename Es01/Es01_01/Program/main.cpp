/********************************************************************************************************************************
 * Numerical Simulation Laboratory - Exercise 01.1
 * Tests of the RANNYU pseudo-random number generator
 *
 * The program generates data for statistical tests on uniform random numbers in [0,1) and writes the results to output files.
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

// Generate random numbers and estimate <r> using the blocking method

	int M=1000000;							// Total number of Monte Carlo throws (number of generated random numbers)
	int N=1000;								// Number of blocks
	int L=M/N;								// Number of throws per block (block length)
	vector<double> rBuffer(L);				// Temporary buffer for values generated in one block
	vector<double> rBlock(N);				// Block averages: rBlock[i] = average of block i (i=0-N-1)
	vector<double> rBlock2(N);				// Squares of block averages: rBlock2[i] = rBlock[i]^2
	vector<double> rMeanValue(N);			// Progressive averages <r>: rMeanValue[i] = <r> estimated over the first (i+1) blocks
	vector<double> rStDevOfTheMean(N);		// Progressive statistical uncertainty of <r>: rStDevOfTheMean[i] = uncertainty of rMeanValue[i]

	for (int i=0; i<N; i++) {
		for (int j=0; j<L; j++) rBuffer[j]=random.Rannyu();			// Generate L random numbers for the i-th block
		rBlock[i]=Average(rBuffer,0,L);								// Compute the average of the i-th block
		rBlock2[i]=rBlock[i]*rBlock[i];
		if (i==0) {
			rMeanValue[i]=rBlock[i];
			rStDevOfTheMean[i]=0.0;									// After the first block the statistical uncertainty is not defined
		} else {
			rMeanValue[i]=Average(rBlock,0,i+1);					// Progressive average over the first (i+1) blocks
			rStDevOfTheMean[i]=sqrt((Average(rBlock2,0,i+1)-rMeanValue[i]*rMeanValue[i])/i);	// Statistical uncertainty of the mean, over the first (i+1) blocks
		}
	}

// Generate random numbers and estimate var=variance using the blocking method - Here var is calculated as <r^2>-<r>^2; the Monte Carlo step consists in squaring the pseudo-random number returned by the RNG

	// Keep same M, N, L
	vector<double> r2Buffer(L);				// Temporary buffer for squared values generated in one block
	vector<double> r2Block(N);				// Block averages: r2Block[i] = average of block i (i=0-N-1)
	vector<double> r2Block2(N);				// Squares of block averages: r2Block2[i] = r2Block[i]^2
	vector<double> r2MeanValue(N);			// Progressive averages <r^2>: r2MeanValue[i] = <r^2> estimated over the first (i+1) blocks
	vector<double> r2StDevOfTheMean(N);		// Progressive statistical uncertainty of <r^2>: r2StDevOfTheMean[i] = uncertainty of r2MeanValue[i]
	vector<double> varMeanValue(N);			// Progressive averages <var>=<r^2>-<r>^2: varMeanValue[i] = r2MeanValue[i]-rMeanValue[N-1]
	vector<double> varStDevOfTheMean(N);	// Progressive statistical uncertainty of <var>: varStDevOfTheMean[i] = uncertainty of varMeanValue[i]

	// Calculate mean and standard deviation of the mean for r^2
	for (int i=0; i<N; i++) {
		for (int j=0; j<L; j++) r2Buffer[j]=pow(random.Rannyu(),2);	// Generate and square L random numbers for the i-th block
		r2Block[i]=Average(r2Buffer,0,L);							// Compute the average of the i-th block: get (r^2)_i
		r2Block2[i]=r2Block[i]*r2Block[i];
		if (i==0) {
			r2MeanValue[i]=r2Block[i];
			r2StDevOfTheMean[i]=0.0;
		} else {
			r2MeanValue[i]=Average(r2Block,0,i+1);					// Progressive average over the first (i+1) blocks: get progressive <r^2>
			r2StDevOfTheMean[i]=sqrt((Average(r2Block2,0,i+1)-r2MeanValue[i]*r2MeanValue[i])/i);	// Statistical uncertainty of the progressive <r^2>, over the first (i+1) blocks
		}
	}

	// Compute progressive averages and statistical uncertainties for var=<r^2>-<r>^2
	varMeanValue[0]=r2MeanValue[0]-rMeanValue[N-1]*rMeanValue[N-1];
	varStDevOfTheMean[0]=0.0;
	for (int i=1; i<N; i++) {
		varMeanValue[i]=r2MeanValue[i]-rMeanValue[N-1]*rMeanValue[N-1];											// Progressive average over the first (i+1) blocks: get progressive <var>
		varStDevOfTheMean[i]=sqrt(pow(r2StDevOfTheMean[i],2)+pow(2.0*rMeanValue[N-1]*rStDevOfTheMean[N-1],2));	// Statistical uncertainty of the progressive <var>, over the first (i+1) blocks - This formula uses the propagation of errors
	}

// Save results to file

	ofstream out("../Output/mean_variance_test.dat");
	if (!out) {
		cerr<<"Output error: Unable to open test.dat"<<endl;
		return 1;
	}

	// Header
	out<<"# Test 1-2: Mean value and variance of uniform random numbers in [0,1)\n";		// Lines starting with '#' are comments and are ignored by numpy.loadtxt in Python
	out<<"# Total number of throws M = "<<M<<"\n";
	out<<"# Number of blocks N = "<<N<<"\n";
	out<<"# Block length L = M/N = "<<L<<"\n";
	out<<"#\n";
	out<<"# Columns:\n";
	out<<"# 1) Block index\n";
	out<<"# 2) Block average\n";
	out<<"# 3) Progressive mean value: <r>_n\n";
	out<<"# 4) Progressive statistical uncertainty of the mean: sigma_<r>_n\n";
	out<<"# 5) Progressive variance: <var>_n\n";
	out<<"# 6) Progressive statistical uncertainty of the variance: sigma_<var>_n\n";
	out<<"# 7) Progressive statistical uncertainty of the mean of the squared values: sigma_<r^2>_n  \n";
	out<<"#\n";

	// Data
	for (int i=0; i<N; i++) out<<i+1<<" "<<rBlock[i]<<" "<<rMeanValue[i]<<" "<<rStDevOfTheMean[i]<<" "<<varMeanValue[i]<<" "<<varStDevOfTheMean[i]<<" "<<r2StDevOfTheMean[i]<<"\n";
	out.close();

// Perform the Pearson's chi^2 test to check the uniformity of the RNG

	int N_tests=1000;										// Number of independent chi^2 tests
	M=10000;												// Number of random numbers generated in each test
	int N_bins=100;											// Number of identical bins dividing the interval [0,1)
	vector<double> chi2(N_tests);							// chi2[j] = chi^2 value of the j-th test
	vector<int> counts(N_bins);								// counts[i] = number of values falling into the i-th bin (for one test)
	double interval_min=0.0;								// Lower bound of the interval
	double interval_max=1.0;								// Upper bound of the interval
	double bin_width=(interval_max-interval_min)/N_bins;	// Width of each bin
	double expectedCountsPerBin=double(M)/double(N_bins);	// Expected number of events per bin for a uniform distribution

	for (int n_test=0; n_test<N_tests; n_test++) {
		chi2[n_test]=0.0;									// Initialize chi^2 for the current test
		fill(counts.begin(),counts.end(),0);				// Reset bin counters
		// Generate M random numbers and fill the histogram
		for (int i=0; i<M; i++)	counts[static_cast<int>(floor((random.Rannyu()-interval_min)/bin_width))]++;
		// Compute chi^2 for the current test
		for (int n_bin=0; n_bin<N_bins; n_bin++) chi2[n_test]+=(counts[n_bin]-expectedCountsPerBin)*(counts[n_bin]-expectedCountsPerBin)/expectedCountsPerBin;
	}
		
// Save chi^2 test results to file

	ofstream out_chi2("../Output/chi2_test.dat");
	if (!out_chi2) {
		cerr<<"Output error: Unable to open chi2_test.dat" << endl;
		return 1;
	}

	// Header
	out_chi2<<"# Pearson chi^2 test for uniform random numbers in [0,1)\n";
	out_chi2<<"# Number of tests = "<<N_tests<<"\n";
	out_chi2<<"# Numbers per test = "<<M<<"\n";
	out_chi2<<"# Number of bins = "<<N_bins<<"\n";
	out_chi2<<"# Degrees of freedom = "<<N_bins-1<<"\n";
	out_chi2<<"#\n";
	out_chi2<<"# Columns:\n";
	out_chi2<<"# 1) Test index\n";
	out_chi2<<"# 2) chi^2 value\n";
	out_chi2<<"#\n";

	// Data
	for (int i=0; i<N_tests; i++) out_chi2<<i+1<<" "<<chi2[i]<<"\n";
	out_chi2.close();

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

