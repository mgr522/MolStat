/*
This work is licensed under the Creative Commons Attribution 3.0 United States
License. To view a copy of this license, visit
http://creativecommons.org/licenses/by/3.0/us/ or send a letter to Creative
Commons, 444 Castro Street, Suite 900, Mountain View, California, 94041, USA.

Copyright (C) 2013 Oak Ridge National Laboratory
*/
/**
 * \file main-simulator.cc
 * \brief Main function for simulating conductance data using Landauer theory.
 *
 * Statistical parameters (for example, the average level energy) are provided
 * on the command-line and used to simulate conductance data. This data can
 * subsequently be binned into a histogram to test the fitting procedures.
 *
 * Two models are presently implemented:
 *    - Symmetric couplings (same coupling used for the left and right leads).
 *    - Asymmetric couplings (different couplings for the two leads).
 *
 * There are seven or eight required command-line arguments, depending on the
 * model:
 *    -# The model to use (`s' for symmetric, `a' for asymmetric).
 *    -# The number of conductance data points to simulate.
 *    -# The Fermi level of the system (eV)
 *    -# The standard deviation in site level energy (eV)
 *    -# The average site level energy (eV)
 *    -# The standard deviation in electrode-channel coupling (eV)
 *    -# The average coupling to one electrode (asymmetric model) or to both
 *       electrodes (symmetric model), in eV.
 *    -# The average coupling to the other electrode (eV). Only required if the
 *       asymmetric model is used.
 *
 * \author Patrick D.\ Williams and Matthew G.\ Reuter
 * \date July 2012, May 2013
 */

#include <cstdio>
#include <cstdlib>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

/**
 * \brief Samples from a normal distribution with the given mean and standard
 *        deviation.
 *
 * \param[in] mean The distribution's average.
 * \param[in] stdev The distribution's standard deviation.
 * \param[in] r The internal workspace for GSL's random number generation.
 * \return The random number.
 */
double normal_random_variable(double mean, double stdev, gsl_rng *r);

/**
 * \brief Evalutates the transmission (Landauer theory) using the specified
 *        physical parameters.
 * 
 * \param[in] gamma1 One channel-lead coupling.
 * \param[in] gamma2 The other channel-lead coupling.
 * \param[in] epsilon The channel's level energy.
 * \param[in] E The Fermi level.
 * \return The transmission.
 */
double transmission(double gamma1, double gamma2, double epsilon, double E);

/**
 * \brief Main function for simulating a histogram.
 *
 * Parses the input parameters and outputs randomly-generated conductance
 * data.
 *
 * \param[in] argc The number of command-line arguments.
 * \param[in] argv The command-line arguments.
 * \return Exit status; 0 for normal.
 */
int main(int argc, char **argv) {
	const gsl_rng_type *T;
	gsl_rng *r;
     
	int i, n;
	double gammaL, gammaR, epsilon, trans;
	double gamma1, gamma2;
	double dgamma;
	double epsilon0;
	double depsilon;
	double E;
	char model;

	if(argc != 8 && argc != 9) {
		fprintf(stderr, "Usage error: ./simulator model n EF deps eps0 dgamma " \
			"gamma1 gamma2 n\n" \
			"   model is the model to use ('s' for symmetric, 'a' for " \
				"asymmetric)\n" \
			"   n is the number of trials\n" \
			"   EF is the Fermi level (eV)\n" \
			"   deps is the standard deviation in site level energy (eV)\n" \
			"   eps0 is the average site level energy (eV)\n" \
			"   dgamma is the standard deviation in the coupling (eV)\n" \
			"   gamma1 is the average coupling for one electrode (eV)\n" \
			"   gamma2 is the average coupling for the other electrode (eV)\n" \
			"\n   NOTE: gamma2 is ignored (and unnecessary) if model == 's'.\n");
		return 0;
	}

	model = argv[1][0];
	if(model != 'a' && model != 's') {
		fprintf(stderr, "Error: Unknown model. Must be 's' or 'a'.\n");
		return 0;
	}

	n = atoi(argv[2]);
	E = atof(argv[3]);
	depsilon = atof(argv[4]);
	epsilon0 = atof(argv[5]);
	dgamma = atof(argv[6]);
	gamma1 = atof(argv[7]);
	if(model == 'a') {
		if(argc == 9)
			gamma2 = atof(argv[8]);
		else {
			fprintf(stderr, "Error: Not enough parameters for the asymmetric " \
				"model.\n       Run ./simulator for the correct usage.\n");
			return 0;
		}
	}

	if(depsilon <= 0.0 || dgamma <= 0.0) {
		fprintf(stderr, "Error: standard deviations must be positive.\n");
		return 0;
	}

	if(gamma1 <= 0.0) {
		fprintf(stderr, "Error: gamma1 must be positive.\n");
		return 0;
	}

	if(model == 'a' && gamma2 <= 0.0) {
		fprintf(stderr, "Error: gamma2 must be positive.\n");
		return 0;
	}

	if(n <= 0) {
		fprintf(stderr, "Error: There must be at least one trial.\n");
		return 0;
	}

	if(gamma1 / dgamma < 4.0 || (model == 'a' && gamma2 / dgamma < 4.0)) {
		fprintf(stderr, "Warning: The models assume gamma1(2) / dgamma >> 0; " \
			"bigger than 4, in practice.\n");
	}

	// Setup the GSL random number generator
	gsl_rng_env_setup();
	T = gsl_rng_default;
	r = gsl_rng_alloc(T);

	// seed the generator
	gsl_rng_set(r, 0xFEEDFACE);

	// Get the requested number of transmission values
	for (i = 0; i < n; ++i) {
		gammaL = normal_random_variable(gamma1, dgamma, r);
		epsilon = normal_random_variable(epsilon0, depsilon, r);

		if(model == 's') {
			trans = transmission(gammaL, gammaL, epsilon, E);
		}
		else if(model == 'a') {
			gammaR = normal_random_variable(gamma2, dgamma, r);
			trans = transmission(gammaL, gammaR, epsilon, E);
		}
		else {
			fprintf(stderr, "An impossible error has occurred.\n");
			trans = 0.0;
		}

		printf("%.6f\n", trans);
	}

	gsl_rng_free(r);
	return 0;
}

double normal_random_variable(double mean, double stdev, gsl_rng *r) {
	return gsl_ran_gaussian(r, stdev) + mean;
}

double transmission(double gamma1, double gamma2, double epsilon, double E) {
	return 4.0*gamma1*gamma2 / 
		(4.0*(E-epsilon)*(E-epsilon) + (gamma1 + gamma2)*(gamma1 + gamma2));
}
