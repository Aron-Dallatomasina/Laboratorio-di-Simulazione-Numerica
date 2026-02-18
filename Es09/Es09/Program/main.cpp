/********************************************************************************************************************************
 * Numerical Simulation Laboratory - Exercise 09.1
 * Genetic Algorithm for the Traveling Salesman Problem (TSP)
 *
 * The program implements a Genetic Algorithm to optimize the path of a salesman visiting N cities exactly once and returning to the starting city.
 * Each chromosome encodes a permutation of cities. The population evolves through rank selection, crossover and mutation operators that preserve the permutation constraint.
 *
 * The algorithm is applied to cities randomly distributed on a circumference or inside a square.
 * The best path length and the average length of the best half population are recorded during evolution.
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
#include "genetic.h"

constexpr double pi=3.14159265358979323846;
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

// Simulation parameters

	// Problem parameters
	int		NumCities;								// Number of cities
	int		NumPop;									// Population size (number of chromosomes)
	int		NumGen;									// Number of generations (GA iterations)
	// Genetic algorithm parameters
	bool	OnCircle;								// true: cities on circumference, false: cities in square
	double	p_sel;									// Rank-selection exponent (selection pressure)
	double	p_cross;								// Crossover probability
	double	p_mut_swap;								// Probability of swap mutation
	double	p_mut_inv;								// Probability of inversion mutation
	double	p_mut_shift;							// Probability of shift mutation
	double	p_mut_mperm;							// Probability of block-permutation mutation

	// Input
	ifstream inChoice("../Input/input.dat");
	if (!inChoice.is_open()) {
		cerr<<"Input error: Unable to open input.dat"<<endl;
		return 1;
	} 	
	inChoice>>NumCities>>NumPop>>NumGen>>OnCircle>>p_sel>>p_cross>>p_mut_swap>>p_mut_inv>>p_mut_shift>>p_mut_mperm; 
	inChoice.close();

	// Generate city positions
	vector<City> Cities;
	if (OnCircle==true)	GenerateCities_Circ(Cities,NumCities,random);				// Generate NumCities cities on a circle
	else				GenerateCities_Square(Cities,NumCities,random);				// Generate NumCities cities in a square
		
	// Initialize the population of chromosomes
	Population pop=InitializePopulation(NumPop,NumCities,Cities,random);
	SortPopulation(pop);															// Sort population by fitness

	// Create/open output files
	ofstream out_length("../Output/best_length.dat");
	if(!out_length.is_open()) {
			cerr<<"Error: unable to open best_length.dat"<<endl;
			return 1;
	}
	ofstream out_half_length("../Output/avg_best_half_length.dat");
	if(!out_half_length.is_open()) {
			cerr<<"Error: unable to open avg_best_half_length.dat"<<endl;
			return 1;
	}
	ofstream out_path("../Output/best_path.dat");
	if(!out_path.is_open()) {
			cerr<<"Error: unable to open best_path.dat"<<endl;
			return 1;
	}
	ofstream out_output("../Output/output.dat");
	if(!out_output.is_open()) {
			cerr<<"Error: unable to open output.dat"<<endl;
			return 1;
	}
	out_length<<"#      	 GEN:      LENGTH:" << endl;
	out_half_length<<"#      	 GEN:      LENGTH:" << endl;
	out_path<<"#      	CITY:			X:     		 Y:" << endl;
	out_output<<"# NumCities = "<<NumCities<<"\n";
	out_output<<"# NumPop = "<<NumPop<<"\n";
	out_output<<"# NumGen = "<<NumGen<<"\n";
	out_output<<"# OnCircle = "<<OnCircle<<"\n";
	out_output<<"# p_sel = "<<p_sel<<"\n";
	out_output<<"# p_cross = "<<p_cross<<"\n";
	out_output<<"# p_mut_swap = "<<p_mut_swap<<"\n";
	out_output<<"# p_mut_inv = "<<p_mut_inv<<"\n";
	out_output<<"# p_mut_shift = "<<p_mut_shift<<"\n";
	out_output<<"# p_mut_mperm = "<<p_mut_mperm<<"\n";
	
	// Loop over NumGen generations
	for (int i=0; i<NumGen; i++) {
		EvolveOneGen(pop,Cities,random,p_sel,p_cross,p_mut_swap,p_mut_inv,p_mut_shift,p_mut_mperm);		// Evolve population by one generation
		out_length<<setw(13)<<(i+1)<<setw(13)<<BestLength(pop)<<endl;									// Best path length at current generation
		out_half_length<<setw(13)<<(i+1)<<setw(13)<<AvgBestHalfLength(pop)<<endl;						// Average path length over best half of the population
		const Chromosome& best=pop.individuals[0];														// Best chromosome (population sorted by fitness)
		for (int j=0; j<NumCities; j++) {
			int idx=best.path[j]-1;																		// Convert 1-based city index to 0-based
			out_path<<setw(13)<<(j+1)<<setw(13)<<Cities[idx].x<<setw(13)<<Cities[idx].y<<endl;			// Cartesian coordinates of ordered best path
		}
	}
	
	out_length.close();
	out_half_length.close();
	out_path.close();
	out_output.close();

// Close the program

	random.SaveSeed();
	return 0; 

}
