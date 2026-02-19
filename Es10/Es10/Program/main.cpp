/********************************************************************************************************************************
 * Numerical Simulation Laboratory - Exercise 10
 * Parallel Genetic Algorithm for the Traveling Salesman Problem (TSP)
 *
 * This program implements a parallel Genetic Algorithm (GA) using MPI to solve the TSP.
 * Each MPI process (a “continent”) performs an independent GA search with its own population and random stream.
 *
 * Every Nmigr generations, the Continents exchange their best individuals, which replace the worst individuals in the receiving populations to enhance global exploration.
 ********************************************************************************************************************************/



// C++ Standard Library
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <iomanip>
#include <cmath>
#include <string>
#include <sstream>
#include <filesystem>
// C++ Standard Template Library
#include <vector>
#include <algorithm>
#include <numeric>
// RNG
#include "random.h"
// GA
#include "genetic.h"
// MPI
#include <mpi.h>

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
	double	p_sel;									// Rank-selection exponent (selection pressure)
	double	p_cross;								// Crossover probability
	double	p_mut_swap;								// Probability of swap mutation
	double	p_mut_inv;								// Probability of inversion mutation
	double	p_mut_shift;							// Probability of shift mutation
	double	p_mut_mperm;							// Probability of block-permutation mutation
	// MPI parameters
	int		NumProc;								// Total number of processes (continents)
	int		Rank;									// Process identifier
	bool	Migration;								// true: migration between continents
	int		NumMigr;								// Number of generations between every migration


	// Input
	ifstream inChoice("../Input/input.dat");
	if (!inChoice.is_open()) {
		cerr<<"Input error: Unable to open input.dat"<<endl;
		return 1;
	} 	
	inChoice>>NumPop>>NumGen>>p_sel>>p_cross>>p_mut_swap>>p_mut_inv>>p_mut_shift>>p_mut_mperm>>Migration>>NumMigr; 
	inChoice.close();

	// MPI initialization
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &NumProc);
    MPI_Comm_rank(MPI_COMM_WORLD, &Rank);

	for (int i=0; i<1000*Rank; i++) random.Rannyu();	// Each process gets a different random sequence

	// Generate city positions	
	vector<City> Cities;
	// Rank 0 reads city coordinates from input file
    if (Rank==0) {
        LoadCities(Cities, "../Input/cap_prov_ita.dat");
		NumCities=(int)Cities.size();
    } else NumCities=0;
	// Rank 0 broadcasts Ncities to all other ranks
    MPI_Bcast(&NumCities,1,MPI_INT,0,MPI_COMM_WORLD);
    if (Rank!=0) Cities.resize(NumCities);
    // Rank 0 broadcasts city coordinates to all other ranks
    vector<double> xs(NumCities);
    vector<double> ys(NumCities);
    if (Rank==0) {
        for (int i=0;i<NumCities;i++) {
            xs[i]=Cities[i].x;
            ys[i]=Cities[i].y;
        }
    }
    MPI_Bcast(xs.data(),NumCities,MPI_DOUBLE,0,MPI_COMM_WORLD);
    MPI_Bcast(ys.data(),NumCities,MPI_DOUBLE,0,MPI_COMM_WORLD);
    for (int i=0;i<NumCities;i++) {
        Cities[i].x=xs[i];
        Cities[i].y=ys[i];
    }

	// Initialize the population of chromosomes
	Population pop=InitializePopulation(NumPop,NumCities,Cities,random);
	SortPopulation(pop);															// Sort population by fitness

	// Vectors used during migration
	vector<int> send_path(NumCities);				// Store the best path to send  
	vector<int> recv_path(NumCities);				// Store the best path received 

	// Create/open output files
	ostringstream rank_num;
	rank_num<<"../Output/r_"<<Rank;
	string dir_path=rank_num.str();
	filesystem::create_directories(dir_path);
	ofstream out_length(dir_path+"/best_length.dat");
	if(!out_length.is_open()) {
			cerr<<"Error: unable to open best_length.dat"<<endl;
			MPI_Abort(MPI_COMM_WORLD, 1);
	}
	ofstream out_half_length(dir_path+"/avg_best_half_length.dat");
	if(!out_half_length.is_open()) {
			cerr<<"Error: unable to open avg_best_half_length.dat"<<endl;
			MPI_Abort(MPI_COMM_WORLD, 1);
	}
	ofstream out_path(dir_path+"/best_path.dat");
	if(!out_path.is_open()) {
			cerr<<"Error: unable to open best_path.dat"<<endl;
			MPI_Abort(MPI_COMM_WORLD, 1);
	}
	
	// Decide whether to perform a migration
	for (int i=0; i<NumGen; i++) {
		bool Migrate=false;
		if (Migration) {
			if (((i+1)%NumMigr==0)&&(NumProc>1)) Migrate=true;
		}	
		// Indipendent evolution
		if (!Migrate) EvolveOneGen(pop,Cities,random,p_sel,p_cross,p_mut_swap,p_mut_inv,p_mut_shift,p_mut_mperm);
		// Random migration between continents
		if (Migrate) {
			// Rank 0 selects two distinct processes for migration: giver and receiver
			int giver=0, receiver=0;
			if (Rank==0) {
				giver=(int)random.Rannyu(0,NumProc);
					do {
					receiver=(int)random.Rannyu(0,NumProc);
				} while (receiver==giver);
			}     
			MPI_Bcast(&giver,1,MPI_INT,0,MPI_COMM_WORLD);
			MPI_Bcast(&receiver,1,MPI_INT,0,MPI_COMM_WORLD);
			// giver sends its best path (firt index)
			if (Rank==giver) {
				send_path=pop.individuals[0].path; 
				MPI_Send(send_path.data(),NumCities,MPI_INT,receiver,0,MPI_COMM_WORLD);
			}
			// receiver replaces its worst path (last index) with the received one
			if (Rank==receiver) {
				MPI_Recv(recv_path.data(),NumCities,MPI_INT,giver,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
				pop.individuals[NumPop-1].path=recv_path; 
				pop.individuals[NumPop-1].check();								// Validate the new path
				CalcFitness(Cities,pop.individuals[NumPop-1]);					// Compute length and fitness of the new path
				SortPopulation(pop);											// Reorder population
			}
		}
		
		// Best path length
		out_length<<setw(13)<<(i+1)<<setw(13)<<BestLength(pop)<<endl;
		// Average path length of best half population
		out_half_length<<setw(13)<<(i+1)<<setw(13)<<AvgBestHalfLength(pop)<<endl;
		// Coordinates of best path in cartesian coordinates
		const Chromosome& best=pop.individuals[0];
		for (int j=0; j<NumCities; j++) {
			int idx=best.path[j]-1;												// Path is 1-based
			out_path<<setw(13)<<(j+1)<< setw(13)<<Cities[idx].x<<setw(13)<<Cities[idx].y<<endl;
		}
	}
	// Find the best global length among all ranks
	struct {double val; int rank;} in_minloc, out_minloc;
	in_minloc.val=BestLength(pop);
	in_minloc.rank=Rank;
	MPI_Reduce(&in_minloc,&out_minloc,1,MPI_DOUBLE_INT,MPI_MINLOC,0,MPI_COMM_WORLD);

	// Rank 0 writes the parallelization results
	if (Rank==0) {
		ofstream out_output("../Output/output.dat");
		if(!out_output.is_open()) cerr<<"PROBLEM: Unable to open output.dat"<<endl;
		
		out_output<<"# PROBLEM PARAMETERS:\n";
		out_output<<"# NumPop = "<<NumPop<<"\n";
		out_output<<"# NumGen = "<<NumGen<<"\n";
		out_output<<"# MPI PARAMETERS:\n";
		out_output<<"# Migration = "<<Migration<<"\n";
		out_output<<"# NumMigr = "<<NumMigr<<"\n";
		out_output<<"# GENETIC PARAMETERS:\n";
		out_output<<"# p_sel = "<<p_sel<<"\n";
		out_output<<"# p_cross = "<<p_cross<<"\n";
		out_output<<"# p_mut_swap = "<<p_mut_swap<<"\n";
		out_output<<"# p_mut_inv = "<<p_mut_inv<<"\n";
		out_output<<"# p_mut_shift = "<<p_mut_shift<<"\n";
		out_output<<"# p_mut_mperm = "<<p_mut_mperm<<"\n";
		out_output<<"# PARALLELIZATION RESULTS:\n";
		out_output<<"# total number of processes = "<<NumProc<<"\n";
		out_output<<"# best global length = "<<out_minloc.val<<" found by rank "<<out_minloc.rank<<"\n";
		out_output.close();
		
		cout<<"Total number of processes = "<<NumProc<<endl;
		cout<<"Best global length = "<<out_minloc.val<<" found by rank "<<out_minloc.rank<<endl;
	}

	out_length.close();
	out_half_length.close();
	out_path.close();


// Close the program
	if (Rank==0) random.SaveSeed();
	MPI_Finalize();
	return 0; 

}
