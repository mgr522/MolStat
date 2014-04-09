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
 * This code produces voltage-dependent conductance histograms. Two types are
 * presently implemented:
 *    - Static conductance (I/V)
 *    - Differential conductance (dI/dV)
 * Symmetric coupling (same coupling to the left and right leads) is assumed,
 * unless otherwise specified.
 *
 * There are eleven required command-line arguments:
 *    -# The type of conductance histogram to simulate:
 *       - `s' for static conductance
 *       - `d' for differential conductance
 *    -# The model to use for generating the transmission function:
 *       - `i' for the voltage-independent model
 *       - `s' for the single-site voltage-dependent model
 *       - `d' for the double-site voltage-dependent model
 *    -# The number of conductance data points to simulate.
 *    -# The Fermi level of the system (eV)
 *    -# The standard deviation in site level energy (eV)
 *    -# The average site level energy (eV)
 *    -# The standard deviation in electrode-channel coupling (eV)
 *    -# The average coupling to both electrodes, (eV)
 *    -# The minimum voltage (V)
 *    -# The maximum voltage (V)
 *    -# The relative voltage drop between the electrodes
 *
 * \author Gaibo Zhang and Matthew G.\ Reuter
 * \date August 2013; April 2014
 */

#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <complex>

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
 * \brief Differential conductance for the voltage-independent model; symmetric
 *        coupling.
 * 
 * \param[in] V The applied voltage.
 * \param[in] gamma The channel-lead coupling.
 * \param[in] epsilon The channel's level energy.
 * \param[in] eta The relative voltage drop between the two electrodes.
 * \param[in] EF The Fermi energy.
 * \return The differential conductance, in units of G0.
 */
double diff_conductance_i(double V, double gamma, double epsilon, double eta,
	double EF);

/**
 * \brief Static conductance for the voltage-independent model; symmetric
 *        coupling.
 * 
 * \param[in] V The applied voltage.
 * \param[in] gamma The channel-lead coupling.
 * \param[in] epsilon The channel's level energy.
 * \param[in] eta The relative voltage drop between the two electrodes.
 * \param[in] EF The Fermi energy.
 * \return The static conductance, in units of G0.
 */
double static_conductance_i(double V, double gamma, double epsilon, double eta,
	double EF);

/**
 * \brief Differential conductance for the single-site voltage-dependent model;
 *        symmetric coupling.
 * 
 * \param[in] V The applied voltage.
 * \param[in] gamma The channel-lead coupling.
 * \param[in] epsilon The channel's level energy.
 * \param[in] EF The Fermi energy.
 * \return The differential conductance, in units of G0.
 */
double diff_conductance_s(double V, double gamma, double epsilon, double eta,
	double EF);

/**
 * \brief Static conductance for the single-site voltage-dependent model;
 *        symmetric coupling.
 * 
 * \param[in] V The applied voltage.
 * \param[in] gamma The channel-lead coupling.
 * \param[in] epsilon The channel's level energy.
 * \param[in] EF The Fermi energy.
 * \return The static conductance, in units of G0.
 */
double static_conductance_s(double V, double gamma, double epsilon, double eta,
	double EF);

/**
 * \brief Differential conductance for the double-site voltage-dependent model;
 *        symmetric coupling.
 * 
 * \param[in] V The applied voltage.
 * \param[in] gamma The channel-lead coupling.
 * \param[in] epsilon The channel's level energy.
 * \param[in] EF The Fermi energy.
 * \return The differential conductance, in units of G0.
 */
double diff_conductance_d(double V, double gamma, double epsilon, double eta,
	double EF);

/**
 * \brief Static conductance for the double-site voltage-dependent model;
 *        symmetric coupling.
 * 
 * \param[in] V The applied voltage.
 * \param[in] gamma The channel-lead coupling.
 * \param[in] epsilon The channel's level energy.
 * \param[in] EF The Fermi energy.
 * \return The static conductance, in units of G0.
 */
double static_conductance_d(double V, double gamma, double epsilon, double eta,
	double EF);

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
	gsl_rng *r;

	char type;
	int i, n;
	double EF;
	double depsilon;
	double epsilon0;
	double dgamma;
	double gamma0;
	double Vmin, Vmax;
	double eta;
	double gamma, epsilon, V, GV;
	double (*cond)(double, double, double, double, double);

	if (argc != 12) {
		fprintf(stderr, "Usage error: ./simulator-v-2d cond model n EF " \
			"depsilon epsilon0 dgamma gamma0 Vmin Vmax eta\n" \
			"   cond is the conductance type: 's' for static or 'd' for " \
			"differential\n" \
			"   model is the model to use: 'i', 's', or 'd'\n" \
			"   n is the number of trials\n" \
			"   EF is the Fermi level (eV)\n" \
			"   depsilon is the standard deviation in site level energy (eV)\n" \
			"   epsilon0 is the average site level energy (eV)\n" \
			"   dgamma is the standard deviation in the coupling (eV)\n" \
			"   gamma0 is the average coupling for one electrode (eV)\n" \
			"   Vmin is the lower bound of the applied bias range (V)\n" \
			"   Vmax is the upper bound of the applied bias range (V)\n" \
			"   eta is the relative voltage drop for one electrode\n" \
			"\n   NOTE: symmetric coupling is assumed.\n");
		return 0;
	}

	// type
	type = *argv[1];
	if(type != 'd' && type != 's') {
		fprintf(stderr, "Error: %c is not a valid conductance type.\n", type);
		return 0;
	}

	// model
	switch(*argv[2]) {
	case 'i':
		if(type == 'd')
			cond = diff_conductance_i;
		else
			cond = static_conductance_i;
		break;
	case 's':
		if(type == 'd')
			cond = diff_conductance_s;
		else
			cond = static_conductance_s;
		break;
	case 'd':
		if(type == 'd')
			cond = diff_conductance_d;
		else
			cond = static_conductance_d;
		break;
	default:
		fprintf(stderr, "Error: Unknown model: '%c'.\n", *argv[1]);
		return 0;
	}
	n = atoi(argv[3]);
	EF = atof(argv[4]);
	depsilon = atof(argv[5]);
	epsilon0 = atof(argv[6]);
	dgamma = atof(argv[7]);
	gamma0 = atof(argv[8]);
	Vmin = atof(argv[9]);
	Vmax = atof(argv[10]);
	eta = atof(argv[11]);

	if (depsilon <= 0.0 || dgamma <= 0.0) {
		fprintf(stderr, "Error: standard deviations must be positive.\n");
		return 0;
	}

	if (gamma0 <= 0.0) {
		fprintf(stderr, "Error: gamma0 must be positive.\n");
		return 0;
	}

	if (n <= 0) {
		fprintf(stderr, "Error: There must be at least one trial.\n");
		return 0;
	}

	if (gamma0 / dgamma < 4.0) {
		fprintf(stderr, "Warning: The model assumes gamma0 / dgamma >> 0; " \
			"bigger than 4, in practice.\n");
	}
	
	if (Vmin > Vmax) {
		fprintf(stderr, "Error: Vmin is lower bound, Vmax is upper bound; " \
			"Vmax > Vmin.\n");
		return 0;
	}

	if (eta > 1.0 || eta < 0.0) {
		fprintf(stderr, "Error: eta is a relative voltage drop on one side; " \
			"0 <= eta <= 1.\n");
		return 0;
	}

	// Setup the GSL random number generator
	gsl_rng_env_setup();
	r = gsl_rng_alloc(gsl_rng_default);

	// Seed the generator
	gsl_rng_set(r, 0xFEEDFACE);

	// Get the requested number of voltage-transmission sets
	for (i = 0; i < n; ++i) {
		V = Vmin + (Vmax - Vmin) * gsl_rng_uniform(r);
		gamma = normal_random_variable(gamma0, dgamma, r);
		epsilon = normal_random_variable(epsilon0, depsilon, r);

		GV = cond(V, gamma, epsilon, eta, EF);

		printf("%.6f %.6f\n", V, GV);
	}

	gsl_rng_free(r);
	return 0;
}

double normal_random_variable(double mean, double stdev, gsl_rng *r) {
	return gsl_ran_gaussian(r, stdev) + mean;
}

// Voltage-independent model
static double transmission_i(double gamma, double epsilon, double E) {
	return gamma*gamma / 
		((E-epsilon)*(E-epsilon) + gamma*gamma);
}

double diff_conductance_i(double V, double gamma, double epsilon, double eta,
	double EF) {

	return eta*transmission_i(gamma, epsilon, EF + eta*V) +
		(1.-eta)*transmission_i(gamma, epsilon, EF + (eta-1.)*V);
}

double static_conductance_i(double V, double gamma, double epsilon, double eta,
	double EF) {

	return gamma/V * (atan((EF-epsilon+eta*V)/gamma)
		- atan((EF-epsilon+(eta-1.)*V)/gamma));
}

// Single-site, voltage-dependent model
static double transmission_s(double V, double gamma, double epsilon, double E)
{
	return gamma*gamma /
		((E-epsilon-V)*(E-epsilon-V) + gamma*gamma);
}

double diff_conductance_s(double V, double gamma, double epsilon, double eta,
	double EF) {

	return (eta-1.)*transmission_s(V, gamma, epsilon, EF + eta*V) +
		(2.-eta)*transmission_s(V, gamma, epsilon, EF + (eta-1.)*V);
}

double static_conductance_s(double V, double gamma, double epsilon, double eta,
	double EF) {

	return gamma/V * (atan((EF-epsilon-V+eta*V)/gamma)
		- atan((EF-epsilon-V+(eta-1.)*V)/gamma));
}

// Double-site, voltage-dependent model
static double transmission_d(double V, double gamma, double epsilon,
	double beta, double E) {

	double temp = 4.*(E-epsilon)*(E-epsilon) - 4.*beta*beta - gamma*gamma - V*V;

	return 16.*gamma*gamma*beta*beta / (temp*temp +
		16.*gamma*gamma*(E-epsilon)*(E-epsilon));
}

static double dtdvint_d(double V, double gamma, double beta, double z) {
	const double bv = 4.*beta*beta + V*V;
	const double bvg = bv + gamma*gamma;
	const std::complex<double> arctan = atan(2.*z /
		std::complex<double>(gamma, -sqrt(bv)));

	return 8.*V*gamma*gamma*beta*beta*z*(4.*z*z + gamma*gamma - 3.*bv) /
		(bv*bvg*(16.*z*z*z*z + 8.*(gamma*gamma - bv)*z*z + bvg*bvg))

		- 8.*V*gamma*beta*beta / (bvg*bvg) * std::real(arctan)

		- 4.*V*gamma*gamma*beta*beta*(3.*bv + gamma*gamma) /
			(bvg*bvg*bv*sqrt(bv)) * std::imag(arctan);
}

double diff_conductance_d(double V, double gamma, double epsilon, double eta,
	double EF) {

	const double beta = -3.0;

	return eta*transmission_d(V, gamma, epsilon, beta, EF + eta*V) +
		(1.-eta)*transmission_d(V, gamma, epsilon, beta, EF + (eta-1.)*V) +
		dtdvint_d(V, gamma, beta, EF - epsilon + eta*V) -
		dtdvint_d(V, gamma, beta, EF - epsilon + (eta-1.)*V);
}

static double tint_d(double V, double gamma, double beta, double z) {
	const std::complex<double> arg(gamma, -sqrt(4.*beta*beta + V*V));

	return 4.*beta*beta*gamma/(V * sqrt(4.*beta*beta+V*V)) *
		std::imag(atan(2.*z/arg)/arg);
}

double static_conductance_d(double V, double gamma, double epsilon, double eta,
	double EF) {

	const double beta = -3.0;

	return tint_d(V, gamma, beta, EF - epsilon + eta*V)
		- tint_d(V, gamma, beta, EF - epsilon + (eta - 1.)*V);
}
