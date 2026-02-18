#pragma once

#include <iostream>
#include <cmath>
#include <algorithm>
#include <vector>
#include "random.h"

using namespace std;


// Structs ------------------------------------------------------------------------------------------------------------------------------------------------------------

// Position of a city (cartesian coordinates)
struct City {
	double x, y;
};

// Chromosome
struct Chromosome {
	vector<int> path;				// Permutation of cities with city 1 fixed at the head
	double fitness=0.0;				// Selection score: 1/L(2)
	double length=0.0;				// Total path length: L(1)
	void check() const;				// Validate the chromosome: path[0]==1 and path is a true permutation
};

// Population of chromosomes
struct Population {
	vector<Chromosome> individuals;
};


// Functions ----------------------------------------------------------------------------------------------------------------------------------------------------------

// City generation
// Generate N cities uniformly on the unit circle
void GenerateCities_Circ(vector<City>& cities, int N, Random& random);
// Generate N cities uniformly in the square
void GenerateCities_Square(vector<City>& cities, int N, Random& random);

// Chromosome
// Create a random chromosome of size N with city 1 fixed at the head
Chromosome BuildRandomChromosome(int N, Random& random);
// Compute length and fitness of a chromosome
void CalcFitness(vector<City>& cities, Chromosome& c);

// Population
// Initialize a population of size M with random chromosomes
Population InitializePopulation(int M, int N, vector<City>& cities, Random& random);
// Compare two chromosomes by fitness
bool CompareChromosomes(const Chromosome& a, const Chromosome& b);
// Sort a population by fitness (best to worst)
void SortPopulation(Population& pop);

// Selection
// Select one index by rank selection
int SelRank(Population& pop, Random& random, double p);
// Select two distinct parent indices by rank selection
void Select2Parents(Population& pop, Random& random, double p, int& i1, int& i2);

// Mutations
// Swap two positions in the tail
void Mutation_Swap(Chromosome& c, Random& random);
// Invert a contiguous segment in the tail
void Mutation_Inversion(Chromosome& c, Random& random);
// Shift a contiguous block in the tail to the right by n positions
void Mutation_Shift(Chromosome& c, Random& random);
// Swap two disjoint blocks of equal length m in the tail
void Mutation_mPermutation(Chromosome& c, Random& random);

// Crossover
// Perform order-based crossover on a segment
void crossover(Chromosome& child, Chromosome& parent2, int pos, int len);

// Evolution, metrics, others
// Evolve the population by one generation
void EvolveOneGen(Population& pop, vector<City>& cities, Random& random,double psel, double pcross, double pm_swap, double pm_inv, double pm_shift, double pm_mperm);
// Get the best path length in the population
double BestLength(Population& pop);
// Get the average path length over the best half of the population
double AvgBestHalfLength(Population& pop);
