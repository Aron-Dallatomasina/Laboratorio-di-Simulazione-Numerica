/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************
 Minor comments added by Aron Dallatomasina, 2026.*/


#pragma once

constexpr double pi=3.14159265358979323846;

// Random numbers
#include "random.h"
Random		random;														// Instance of the random number generator
int			seed[4];													// Seed

// Parameters, observables
const int	m_props=1000;												// Maximum number of measurable observables (upper bound used to size arrays for block and global averages)
int			n_props;													// Actual number of observables used in the simulation
int			iu, ic, im, ix, ig;											// Indices of the observables in the walker array (energy, heat capacity, magnetization, magnetic susceptibility, radial distribution function g(r))
double		nbins;														// Number of bins used for histogram-based observables (e.g. g(r)); not used in the present Ising 1D simulation
double		walker[m_props];											// Instantaneous values of the observables at each simulation step

// Averages
double		blk_av[m_props], blk_norm, accepted, attempted;				// Block averages, normalization factor, Metropolis statistics
double		glob_av[m_props], glob_av2[m_props];						// Progressive global averages and squared averages (for error evaluation)
double		stima_u, stima_c, stima_m, stima_x, stima_g;				// Block estimates of physical observables
double		err_u, err_c, err_m, err_x, err_g;							// Statistical uncertainties (data blocking)

// Configuration
const int	m_spin=50;													// Maximum number of spins allowed in the simulation
double		s[m_spin];													// Current spin configuration

// Thermodynamical state
int			nspin;														// Number of spins in the system
double		beta, temp, J, h;											// Thermodynamic parameters (1/T, T, spin coupling constant, external magnetic field)

// Simulation
int			nstep, nblk, metro, restart;								// Steps per block, number of blocks, Monte Carlo algorithm selector (Gibbs/Metropolis), restart flag

// Functions
void		Input(void);												// Read input parameters and initialize the system
void		Reset(int);													// Reset block statistics
void		Accumulate(void);											// Accumulate observables within a block
void		Averages(int);												// Compute block averages and statistical uncertainties
void		Move();														// Perform one Metropolis or Gibbs move
void		ConfFinal(void);											// Save final configuration
void		Measure(void);												// Measure instantaneous physical observables
double		Boltzmann(int,int);											// Compute energy contribution of a trial spin value (Ising model)
int			Pbc(int);													// Apply periodic boundary conditions
double		Error(double,double,int);									// Statistical uncertainty (data blocking)
