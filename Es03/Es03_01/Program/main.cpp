/********************************************************************************************************************************
 * Numerical Simulation Laboratory - Exercise 03.1
 * Monte Carlo pricing of European options
 *
 * The program estimates the price of European Call and Put options within the Black–Scholes framework by simulating the stochastic evolution of the underlying asset price S(t), modeled as a geometric Brownian motion.
 * 
 * Two different approaches are considered:
 * 	1. Direct sampling of the final asset price S(T)
 * 	2. Discretized simulation of the stochastic path in time
 * Statistical uncertainties are estimated using the data blocking technique.
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
// Simulate the asset price S(t) at time t according to the geometric Brownian motion starting from S0 at time t0, with risk-free rate r and volatility sigma
double S(Random &random, double t, double t0, double S0, double r, double sigma);
// Return the maximum between two real numbers a and b
double Max(double a, double b);

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

// 1. Monte Carlo pricing of European Call and Put options by direct sampling of S(T):
// the asset price S(t) follows a geometric Brownian motion and S(T) is sampled directly. The discounted payoff is computed for each throw and averaged using the data blocking technique.

	// Simulation parameters
	int M=1000000;							// Total number of throws (simulated price evolutions)
	int N=1000;								// Number of blocks. Each block yields an estimate of the option price
	int L=M/N;								// Number of throws per block (block length)
	double t0=0.0;							// Initial time
	double T=1.0;							// Final time (expiry date of the option)
	double S0=100.0;						// Initial spot price of the underlying: S0=S(t0)
	double K=100.0;							// Strike price of the underlying
	double r=0.1;							// Risk-free interest rate
	double sigma=0.25;						// Volatility of the underlying
	double Discount=exp(-r*T);				// Discount factor exp(-rT)

	// Arrays for data storage and blocking
	double ST=0.0;								// Temporary variable for the simulated spot price at time T at a given Monte Carlo throw
	// Call option
	vector<double> callpriceBuffer(L);			// Temporary buffer for call option prices generated in one block
	vector<double> callpriceBlock(N);			// Block averages for call option prices: callpriceBlock[i] = average of block i (i=0-N-1)
	vector<double> callpriceBlock2(N);			// Squares of block averages for call option prices: callpriceBlock2[i] = callpriceBlock[i]^2
	vector<double> callpriceMeanValue(N);		// Progressive averages <callprice>: callpriceMeanValue[i] = <callprice> estimated over the first (i+1) blocks
	vector<double> callpriceStDevOfTheMean(N);	// Progressive statistical uncertainty of <callprice>: callpriceStDevOfTheMean[i] = uncertainty of callpriceMeanValue[i]
	// Put option
	vector<double> putpriceBuffer(L);			// Temporary buffer for put option prices generated in one block
	vector<double> putpriceBlock(N);			// Block averages for put option prices: putpriceBlock[i] = average of block i (i=0-N-1)
	vector<double> putpriceBlock2(N);			// Squares of block averages for put option prices: putpriceBlock2[i] = putpriceBlock[i]^2
	vector<double> putpriceMeanValue(N);		// Progressive averages <putprice>: putpriceMeanValue[i] = <putprice> estimated over the first (i+1) blocks
	vector<double> putpriceStDevOfTheMean(N);	// Progressive statistical uncertainty of <putprice>: putpriceStDevOfTheMean[i] = uncertainty of putpriceMeanValue[i]

	// Simulation
	for (int i=0; i<N; i++) {															// Loop over blocks

		for (int j=0; j<L; j++) {														// Loop over spot price evolutions in a block
			ST=S(random,T,t0,S0,r,sigma);												// Compute (simulated) spot price at time T
			callpriceBuffer[j]=Discount*Max(0.0,ST-K);									// Compute call option price as the discounted call payoff
			putpriceBuffer[j]=Discount*Max(0.0,K-ST);									// Compute put option price as the discounted put payoff
		}

		callpriceBlock[i]=Average(callpriceBuffer,0,L);									// Compute the average of the i-th block
		callpriceBlock2[i]=callpriceBlock[i]*callpriceBlock[i];							// Compute its squared value - used in computing statistical uncertainty
		putpriceBlock[i]=Average(putpriceBuffer,0,L);
		putpriceBlock2[i]=putpriceBlock[i]*putpriceBlock[i];

		if (i==0) {
			callpriceMeanValue[i]=callpriceBlock[i];
			callpriceStDevOfTheMean[i]=0.0;												// After the first block the statistical uncertainty is not defined
			putpriceMeanValue[i]=putpriceBlock[i];
			putpriceStDevOfTheMean[i]=0.0;
		} else {
			callpriceMeanValue[i]=Average(callpriceBlock,0,i+1);						// Progressive average over the first (i+1) blocks
			callpriceStDevOfTheMean[i]=sqrt((Average(callpriceBlock2,0,i+1)-callpriceMeanValue[i]*callpriceMeanValue[i])/i);	// Statistical uncertainty of the mean, over the first (i+1) blocks
			putpriceMeanValue[i]=Average(putpriceBlock,0,i+1);
			putpriceStDevOfTheMean[i]=sqrt((Average(putpriceBlock2,0,i+1)-putpriceMeanValue[i]*putpriceMeanValue[i])/i);
		}

	}

// Save results to file

	ofstream out_direct("../Output/results_direct.dat");
	if (!out_direct.is_open()) {
		cerr<<"Error: unable to open results_direct.dat"<<endl;
		return 1;
	}
	out_direct<<"# European option pricing by Monte Carlo direct sampling of S(T)\n";
	out_direct<<"# Total number of throws (option price simulations): M = "<<M<<"\n";
	out_direct<<"# Number of blocks: N = "<<N<<"\n";
	out_direct<<"# Throws per block: L = "<<L<<"\n";
	out_direct<<"# Initial time: t0 = "<<t0<<"\n";
	out_direct<<"# Final time (expiry date): T = "<<T<<"\n";
	out_direct<<"# Initial spot price of the underlying: S0 = "<<S0<<"\n";
	out_direct<<"# Strike price: K = "<<K<<"\n";
	out_direct<<"# Risk-free interest rate: r = "<<r<<"\n";
	out_direct<<"# Market volatility: sigma = "<<sigma<<"\n";
	out_direct<<"# Discount factor: exp(-rT) = "<<Discount<<"\n";
	out_direct<<"# Columns:\n";
	out_direct<<"# 1) block index\n";
	out_direct<<"# 2) progressive mean of the call option price\n";
	out_direct<<"# 3) progressive statistical uncertainty\n";
	out_direct<<"# 4) progressive mean of the put option price\n";
	out_direct<<"# 5) progressive statistical uncertainty\n";

	for (int i=0; i<N; i++) out_direct<<i+1<<" "<<callpriceMeanValue[i]<<" "<<callpriceStDevOfTheMean[i]<<" "<<putpriceMeanValue[i]<<" "<<putpriceStDevOfTheMean[i]<<"\n";
	out_direct.close();

// 2. Monte Carlo pricing of European Call and Put options by discretized simulation of the spot price path:
// the asset price S(t) follows a geometric Brownian motion and its evolution is simulated through a finite number of time steps up to T. The discounted payoff is computed for each throw and averaged using the data blocking technique.

	// Simulation parameters
	M=1000000;								// Total number of throws (simulated price paths)
	N=1000;									// Number of blocks. Each block yields an estimate of the option price
	L=M/N;									// Number of throws per block (block length)
	t0=0.0;									// Initial time
	T=1.0;									// Final time (expiry date of the option)
	int N_steps=100;						// Number of time steps used to discretize the interval [t0,T]
	double T_step=(T-t0)/double(N_steps);	// Time step size	
	S0=100.0;								// Initial spot price of the underlying: S0=S(t0)
	K=100.0;								// Strike price of the underlying
	r=0.1;									// Risk-free interest rate
	sigma=0.25;								// Volatility of the underlying
	Discount=exp(-r*T);						// Discount factor exp(-rT)

	// Arrays for data storage and blocking
	double St=0.0;							// Temporary variable for the simulated spot price at time t
	// Call option
	callpriceBuffer.assign(L,0.0);			// Temporary buffer for call option prices generated in one block
	callpriceBlock.assign(N,0.0);			// Block averages for call option prices: callpriceBlock[i] = average of block i (i=0-N-1)
	callpriceBlock2.assign(N,0.0);			// Squares of block averages for call option prices: callpriceBlock2[i] = callpriceBlock[i]^2
	callpriceMeanValue.assign(N,0.0);		// Progressive averages <callprice>: callpriceMeanValue[i] = <callprice> estimated over the first (i+1) blocks
	callpriceStDevOfTheMean.assign(N,0.0);	// Progressive statistical uncertainty of <callprice>: callpriceStDevOfTheMean[i] = uncertainty of callpriceMeanValue[i]
	// Put option
	putpriceBuffer.assign(L,0.0);			// Temporary buffer for put option prices generated in one block
	putpriceBlock.assign(N,0.0);			// Block averages for put option prices: putpriceBlock[i] = average of block i (i=0-N-1)
	putpriceBlock2.assign(N,0.0);			// Squares of block averages for put option prices: putpriceBlock2[i] = putpriceBlock[i]^2
	putpriceMeanValue.assign(N,0.0);		// Progressive averages <putprice>: putpriceMeanValue[i] = <putprice> estimated over the first (i+1) blocks
	putpriceStDevOfTheMean.assign(N,0.0);	// Progressive statistical uncertainty of <putprice>: putpriceStDevOfTheMean[i] = uncertainty of putpriceMeanValue[i]

	// Simulation
	for (int i=0; i<N; i++) {															// Loop over blocks

		for (int j=0; j<L; j++) {														// Loop over spot price evolutions in a block
			St=S0;																		// Set St as the initial spot price
			for (int k=0; k<N_steps; k++) St=S(random,T_step,0.0,St,r,sigma);			// Loop over the time steps: simulate the evolution of S(t) through N_steps increments of size T_step
			callpriceBuffer[j]=Discount*Max(0.0,St-K);									// Compute call option price as the discounted call payoff
			putpriceBuffer[j]=Discount*Max(0.0,K-St);									// Compute put option price as the discounted put payoff
		}

		callpriceBlock[i]=Average(callpriceBuffer,0,L);									// Compute the average of the i-th block
		callpriceBlock2[i]=callpriceBlock[i]*callpriceBlock[i];							// Compute its squared value - used in computing statistical uncertainty
		putpriceBlock[i]=Average(putpriceBuffer,0,L);
		putpriceBlock2[i]=putpriceBlock[i]*putpriceBlock[i];

		if (i==0) {
			callpriceMeanValue[i]=callpriceBlock[i];
			callpriceStDevOfTheMean[i]=0.0;												// After the first block the statistical uncertainty is not defined
			putpriceMeanValue[i]=putpriceBlock[i];
			putpriceStDevOfTheMean[i]=0.0;
		} else {
			callpriceMeanValue[i]=Average(callpriceBlock,0,i+1);						// Progressive average over the first (i+1) blocks
			callpriceStDevOfTheMean[i]=sqrt((Average(callpriceBlock2,0,i+1)-callpriceMeanValue[i]*callpriceMeanValue[i])/i);	// Statistical uncertainty of the mean, over the first (i+1) blocks
			putpriceMeanValue[i]=Average(putpriceBlock,0,i+1);
			putpriceStDevOfTheMean[i]=sqrt((Average(putpriceBlock2,0,i+1)-putpriceMeanValue[i]*putpriceMeanValue[i])/i);
		}

	}

// Save results to file

	ofstream out_discrete("../Output/results_discrete.dat");
	if (!out_discrete.is_open()) {
		cerr<<"Error: unable to open results_discrete.dat"<<endl;
		return 1;
	}
	out_discrete<<"# European option pricing by Monte Carlo discretized simulation of the underlying path\n";
	out_discrete<<"# Total number of throws (option price simulations): M = "<<M<<"\n";
	out_discrete<<"# Number of blocks: N = "<<N<<"\n";
	out_discrete<<"# Throws per block: L = "<<L<<"\n";
	out_discrete<<"# Initial time: t0 = "<<t0<<"\n";
	out_discrete<<"# Final time (expiry date): T = "<<T<<"\n";
	out_discrete<<"# Number of time steps: N_steps = "<<N_steps<<"\n";
	out_discrete<<"# Time step size: T_step = "<<T_step<<"\n";
	out_discrete<<"# Initial spot price of the underlying: S0 = "<<S0<<"\n";
	out_discrete<<"# Strike price: K = "<<K<<"\n";
	out_discrete<<"# Risk-free interest rate: r = "<<r<<"\n";
	out_discrete<<"# Market volatility: sigma = "<<sigma<<"\n";
	out_discrete<<"# Discount factor: exp(-rT) = "<<Discount<<"\n";
	out_discrete<<"# Columns:\n";
	out_discrete<<"# 1) block index\n";
	out_discrete<<"# 2) progressive mean of the call option price\n";
	out_discrete<<"# 3) progressive statistical uncertainty\n";
	out_discrete<<"# 4) progressive mean of the put option price\n";
	out_discrete<<"# 5) progressive statistical uncertainty\n";

	for (int i=0; i<N; i++) out_discrete<<i+1<<" "<<callpriceMeanValue[i]<<" "<<callpriceStDevOfTheMean[i]<<" "<<putpriceMeanValue[i]<<" "<<putpriceStDevOfTheMean[i]<<"\n";
	out_discrete.close();


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

double S(Random &random, double t, double t0, double S0, double r, double sigma) {
	double z=random.Gauss(0.0,1.0);
	return S0*exp( (r-0.5*sigma*sigma)*(t-t0) + sigma*z*sqrt(t-t0) );
}

double Max(double a, double b) {
	return (a>b)? a:b;
}
