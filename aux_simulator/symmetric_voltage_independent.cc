/**
 * \file voltage_independent.cc
 * \brief Implementation of the voltage-independent tight-binding model for
 *        calculating conductances.
 *
 * \author Matthew G.\ Reuter
 * \date April 2014
 */

#include "symmetric_voltage_independent.h"
#include <cmath>
#include <string>
#include <vector>
#include <stdexcept>
#include "../string_tools.h"

using namespace std;

shared_ptr<ConductanceModel> SymmetricVoltageIndependentModel::create_model(
	FILE *f) {

	shared_ptr<RandomDistribution> eta, eps, gamma;
	string line;
	vector<string> tokens;
	shared_ptr<RandomDistribution> *dist;

	eta = eps = gamma = nullptr;

	// need to read in three random variables: eta, epsilon, gamma.
	// read in three lines: the first token tells us the variable
	for(int i = 0; i < 3; ++i) {
		// this call might throw a runtime_error, but that's okay.
		line = getline(f);

		tokenize(line, tokens);
		if(tokens.size() == 0)
			throw runtime_error("Blank line encountered.");
		make_lower(tokens[0]);

		// figure out which distribution we're reading
		if(tokens[0] == "eta") {
			printf("Reading eta\n");
			dist = &eta;
		}
		else if(tokens[0] == "epsilon" || tokens[0] == "eps") {
			printf("Reading epsilon\n");
			dist = &eps;
		}
		else if(tokens[0] == "gamma") {
			printf("Reading gamma\n");
			dist = &gamma;
		}
		else
			throw runtime_error("Unrecognized parameter.");

		tokens.erase(tokens.begin());

		try {
			*dist = distribution_from_tokens(tokens);
		}
		catch(const invalid_argument &e) {
			throw runtime_error("Unable to form a random number distribution.");
		}
	}

	if(eta == nullptr || eps == nullptr || gamma == nullptr)
		throw runtime_error("Distributions for eta, epsilon, and/or gamma are " \
			"unspecified.");

	return make_shared<SymmetricVoltageIndependentModel>(eta, eps, gamma);
}

SymmetricVoltageIndependentModel::SymmetricVoltageIndependentModel(
	shared_ptr<const RandomDistribution> eta,
	shared_ptr<const RandomDistribution> eps,
	shared_ptr<const RandomDistribution> gamma)
	: ConductanceModel(eta), dist_eps(eps), dist_gamma(gamma) {}

double SymmetricVoltageIndependentModel::static_conductance(
	shared_ptr<gsl_rng> r, const double EF, const double V) const {

	// get model parameters from the random distributions
	double eta = dist_eta->sample(r);
	double eps = dist_eps->sample(r);
	double gamma = dist_gamma->sample(r);

	return static_conductance(EF, V, eta, eps, gamma);
}

double SymmetricVoltageIndependentModel::diff_conductance(
	shared_ptr<gsl_rng> r, const double EF, const double V) const {

	// get model parameters from the random distributions
	double eta = dist_eta->sample(r);
	double eps = dist_eps->sample(r);
	double gamma = dist_gamma->sample(r);

	return diff_conductance(EF, V, eta, eps, gamma);
}

double SymmetricVoltageIndependentModel::transmission(const double E,
	const double eps, const double gamma) {

	return gamma*gamma / ((E-eps)*(E-eps) + gamma*gamma);
}

double SymmetricVoltageIndependentModel::static_conductance(const double EF,
	const double V, const double eta, const double eps, const double gamma) {

	return gamma / V
		* (atan((EF-eps+eta*V) / gamma) - atan((EF-eps+(eta-1.)*V) / gamma));
}

double SymmetricVoltageIndependentModel::diff_conductance(const double EF,
	const double V, const double eta, const double eps, const double gamma) {

	return eta * transmission(EF + eta*V, eps, gamma) +
		(1.-eta) * transmission(EF + (eta-1.)*V, eps, gamma);
}
