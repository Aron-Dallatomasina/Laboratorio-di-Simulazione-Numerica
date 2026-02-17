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
#include <vector>

constexpr double pi=3.14159265358979323846;

// Random numbers
#include "random.h"
Random		random;														// Instance of the random number generator
int			seed[4];													// Seed

// Parameters, observables
const int	m_props=1000;												// Maximum number of measurable observables (upper bound used to size arrays for block and global averages)
int			n_props;													// Actual number of observables used in the simulation
int			iv, ik, it, ie, ip, ig;										// Indices of the observables in the walker array (potential, kinetic, temperature, total energy, pressure, starting index of g(r))
double		vtail, ptail, bin_size, nbins, sd;							// Variables for tail corrections and radial distribution function g(r) (not used in this version)
double		walker[m_props];											// Instantaneous values of the observables at each simulation step

// Averages
double				blk_av[m_props], blk_norm, accepted, attempted;				// Block averages, normalization factor, Metropolis statistics
double				glob_av[m_props], glob_av2[m_props];						// Progressive global averages and squared averages (for error evaluation)
double				stima_pot, stima_pres, stima_kin, stima_etot, stima_temp;	// Block estimates of physical observables
double				err_pot, err_press, err_kin, err_etot, err_temp;			// Statistical uncertainties (data blocking)
std::vector<double>	stima_g, err_g;												// Block estimates and statistical uncertainties for each bin of the g(r) histogram
double				r_lower, r_upper, deltaV;									// Lower and upper bounds of the radial bin and corresponding shell volume for g(r) normalization


// Configuration
const int	m_part=108;													// Maximum number of particles allowed in the simulation
double		x[m_part], y[m_part], z[m_part];							// Current particle positions
double		xold[m_part], yold[m_part], zold[m_part];					// Previous particle positions (needed for Verlet integration)
double		vx[m_part], vy[m_part], vz[m_part];							// Particle velocities

// Thermodynamical state
int			npart;														// Number of particles in the system
double		beta, temp, energy, vol, rho, box, rcut;					// Thermodynamic parameters (1/T, T, volume, density, box length, cutoff radius)

// Simulation
int			iNVET, nstep, nblk, restart;								// Ensemble selector (NVE/NVT), steps per block, number of blocks, restart flag
double		delta;														// Time step (MD) or maximum displacement (MC)

// Functions
void		Input(void);												// Read input parameters and initialize the system
void		Reset(int);													// Reset block statistics
void		Accumulate(void);											// Accumulate observables within a block
void		Averages(int);												// Compute block averages and statistical uncertainties
void		Move(void);													// Perform one MC or MD move
void		ConfFinal(void);											// Save final configuration and velocities
void		ConfXYZ(int);												// Write configuration in XYZ format
void		Measure(void);												// Measure instantaneous physical observables
double		Boltzmann(double,double,double,int);						// Compute single-particle potential energy (MC)
double		Pbc(double);												// Apply periodic boundary conditions
double		Error(double,double,int);									// Statistical uncertainty (data blocking)
double		Force(int,int);												// Compute force component acting on a particle
