/**
 * \file main-simulator-v-2d.cc
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
 * \author Matthew G.\ Reuter
 * \date April 2014
 */

#include <memory>
#include <cstdio>
#include <stdexcept>
#include <string>
#include <vector>
#include <functional>
#include "string_tools.h"
#include "aux_simulator/rng.h"
#include "aux_simulator/symmetric_voltage_independent.h"

using namespace std;

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

	// initialize the GSL random number generator
	gsl_rng_env_setup();
	shared_ptr<gsl_rng> r(gsl_rng_alloc(gsl_rng_default), &gsl_rng_free);
	gsl_rng_set(r.get(), 0xFEEDFACE);

	char type;
	int i, n;
	double EF;
	double V, GV;
	shared_ptr<ConductanceModel> model;
	shared_ptr<RandomDistribution> dist_V;
	function<shared_ptr<ConductanceModel>(FILE *)> creator;
	function<double(const double, const double)> conductance_function;

	string line;
	vector<string> tokens;

	// setup the simulation -- read in parameters from stdin
	// Line 1: One token specifying the model to use
	try {
		line = getline(stdin);
	}
	catch(const runtime_error &e) {
		fprintf(stderr, "Error: %s.\n", e.what());
		return 0;
	}

	tokenize(line, tokens);
	if(tokens.size() < 1) {
		fprintf(stderr, "Error: model name expected in line 1.\n");
		return 0;
	}

	// make a link to the correct model creator (need to continue reading data,
	// store for later)
	make_lower(tokens[0]);
	if(tokens[0] == "symmetricvoltageindependentmodel")
		creator = SymmetricVoltageIndependentModel::create_model;
	else {
		fprintf(stderr, "Error: Unrecognized model in line 1. Options are:\n" \
			"   VoltageIndependentModel - Voltage-Independent Transmission\n");
		return 0;
	}

	// Line 2: Static or differential conductance?
	try {
		line = getline(stdin);
	}
	catch(const runtime_error &e) {
		fprintf(stderr, "Error: %s.\n", e.what());
		return 0;
	}
	tokenize(line, tokens);
	if(tokens.size() < 1) {
		fprintf(stderr, "Error: conductance type expected in line 2.\n");
		return 0;
	}
	make_lower(tokens[0]);
	if(tokens[0] == "static")
		type = 's';
	else if(tokens[0] == "differential")
		type = 'd';
	else {
		fprintf(stderr, "Error: Unrecognized conductance type in line 2.\n   " \
			"It must be \"Static\" or \"Differential\".\n");
		return 0;
	}

	// Line 3: number of trials
	try {
		line = getline(stdin);
	}
	catch(const runtime_error &e) {
		fprintf(stderr, "Error: %s.\n", e.what());
		return 0;
	}
	tokenize(line, tokens);
	if(tokens.size() < 1) {
		fprintf(stderr, "Error: number of trials expected in line 3.\n");
		return 0;
	}
	try {
		n = stoi(tokens[0]);
	}
	catch(const invalid_argument &e) {
		fprintf(stderr, "Error: unrecognizable number '%s'.\n",
			tokens[0].c_str());
		return 0;
	}

	// Line 4: Fermi level
	try {
		line = getline(stdin);
	}
	catch(const runtime_error &e) {
		fprintf(stderr, "Error: %s.\n", e.what());
		return 0;
	}
	tokenize(line, tokens);
	if(tokens.size() < 1) {
		fprintf(stderr, "Error Fermi energy expected in line 4.\n");
		return 0;
	}
	try {
		EF = stod(tokens[0]);
	}
	catch(const invalid_argument &e) {
		fprintf(stderr, "Error: unable to parse %s to the Fermi energy.\n",
			tokens[0].c_str());
		return 0;
	}

	// Line 5: Random distribution for voltages
	try {
		line = getline(stdin);
	}
	catch(const runtime_error &e) {
		fprintf(stderr, "Error: %s.\n", e.what());
		return 0;
	}
	tokenize(line, tokens);
	if(tokens.size() == 0) {
		fprintf(stderr, "Error: expecting voltage distribution on line 5.\n");
		return 0;
	}
	tokens.erase(tokens.begin());
	try {
		dist_V = distribution_from_tokens(tokens);
	}
	catch(const invalid_argument &e) {
		fprintf(stderr, "Error: unable to form a random number distribution " \
			"from:\n   %s\n%s\n", line.c_str(), e.what());
		return 0;
	}

	// invoke the creator for the particular model to read in any other
	// details
	try {
		model = creator(stdin);
	}
	catch(const runtime_error &e) {
		fprintf(stderr, "Error: invalid model parameters.\n   %s\n", e.what());
		return 0;
	}

	// set the conductance function
	if(type == 's')
		conductance_function = [=] (const double a, const double b) {
			return model->static_conductance(r, a, b);
		};
	else if(type == 'd')
		conductance_function = [=] (const double a, const double b) {
			return model->diff_conductance(r, a, b);
		};

	// Get the requested number of voltage-conductance data points
	for (i = 0; i < n; ++i) {
		V = dist_V->sample(r);
		GV = conductance_function(EF, V);

		printf("%.6f %.6f\n", V, GV);
	}

	return 0;
}

#if 0
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
#endif
