/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************

 * Implementation of the RANNYU random number generator. */


#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include "random.h"
constexpr double pi=3.14159265358979323846;
using namespace std;


// Initialize the generator state and select the random stream
void Random::SetRandom(int *s, int p1, int p2) {
	// Multiplier (a) written in base 4096 (fixed for the RANNYU generator)
	m1=502;
	m2=1521;
	m3=4071;
	m4=2107;
	// Initialize the 48-bit internal state (x_{n}) from the input seed (state written in base 4096: 4 × 12-bit blocks)
	l1=s[0];
	l2=s[1];
	l3=s[2];
	l4=s[3];
	// Addend (c) written in base 4096. Only the least significant 24 bits (n3, n4) are used to define statistically independent random streams
	n1=0;
	n2=0;
	n3=p1;
	n4=p2;
}

// Save the current generator state to file
void Random::SaveSeed() {
	ofstream out("../Output/seed.out");
	if (!out) {
		cerr<<"Output error: Unable to open seed.out"<<endl;
		return;
	}
	out<<"RANDOMSEED\t"<<l1 <<" "<<l2<<" "<<l3<<" "<<l4<<endl;		// The "RANDOMSEED" keyword makes seed.out compatible with seed.in, allowing the simulation to be restarted from the saved state.
}

// Generate a random number uniformly distributed in [0,1). This function advances the 48-bit LCG state using a base-4096 (4 × 12-bit blocks) representation with explicit carry handling
double Random::Rannyu() {
	// Definitions
	const double twom12=0.000244140625;								// 1/4096=2^-12, used to reconstruct the double in [0,1)
	int i1, i2, i3, i4;
	double r;
	// Compute the product a*x_{n}+c in base 4096 (partial products plus addend coefficients)
	i1=l1*m4+l2*m3+l3*m2+l4*m1+n1;
	i2=l2*m4+l3*m3+l4*m2+n2;
	i3=l3*m4+l4*m3+n3;
	i4=l4*m4+n4;
    // Propagate carries and update the 4 × 12-bit state blocks
	l4=i4%4096;
	i3=i3+i4/4096;
	l3=i3%4096;
	i2=i2+i3/4096;
	l2=i2%4096;
	l1=(i1+i2/4096)%4096;
    // Reconstruct the random number in [0,1). r=(l1*4096^3+l2*4096^2+l3*4096+l4)/2^48
    r=twom12*(l1+twom12*(l2+twom12*(l3+twom12*l4)));
    return r;
}

// Generate a random number uniformly distributed in [min, max)
double Random::Rannyu(double min, double max) {
	return min+(max-min)*Rannyu();
}

// Generate a random number with Gaussian distribution N(mean, sigma) using the Box–Muller transform
double Random::Gauss(double mean, double sigma) {
	double s=Rannyu();
	double t=Rannyu();
	double x=sqrt(-2.0*log(1.0-s))*cos(2.0*pi*t);
	return mean+x*sigma;
}

// Generate a random number with exponential distribution with parameter lambda
double Random::Exponential(double lambda) {
	return -log(1-Rannyu())/lambda;
}

// Generate a random number with Cauchy-Lorentz distribution with center mu and width Gamma
double Random::CauchyLorentz(double mu, double Gamma) {
	return mu+Gamma*tan(pi*(Rannyu()-0.5));
}
