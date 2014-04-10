/**
 * \file voltage_independent.cc
 * \brief Implementation of the voltage-independent tight-binding model for
 *        calculating conductances.
 *
 * \author Matthew G.\ Reuter
 * \date April 2014
 */

#include "voltage_independent.h"
#include <cmath>

VoltageIndependentModel::VoltageIndependentModel(
	shared_ptr<const RandomDistribution> eta,
	shared_ptr<const RandomDistribution> eps,
	shared_ptr<const RandomDistribution> gamma)
	: ConductanceModel(eta), dist_eps(eps), dist_gamma(gamma) {}

double VoltageIndependentModel::static_conductance(shared_ptr<gsl_rng> r,
	const double EF, const double V) const {

	// get model parameters from the random distributions
	double eta = dist_eta->sample(r);
	double eps = dist_eps->sample(r);
	double gamma = dist_gamma->sample(r);

	return static_conductance(EF, V, eta, eps, gamma);
}

double VoltageIndependentModel::diff_conductance(shared_ptr<gsl_rng> r,
	const double EF, const double V) const {

	// get model parameters from the random distributions
	double eta = dist_eta->sample(r);
	double eps = dist_eps->sample(r);
	double gamma = dist_gamma->sample(r);

	return diff_conductance(EF, V, eta, eps, gamma);
}

double VoltageIndependentModel::transmission(const double E, const double eps,
	const double gamma) {

	return gamma*gamma / ((E-eps)*(E-eps) + gamma*gamma);
}

double VoltageIndependentModel::static_conductance(const double EF,
	const double V, const double eta, const double eps, const double gamma) {

	return gamma / V
		* (atan((EF-eps+eta*V) / gamma) - atan((EF-eps+(eta-1.)*V)));
}

double VoltageIndependentModel::diff_conductance(const double EF,
	const double V, const double eta, const double eps, const double gamma) {

	return eta * transmission(EF + eta*V, eps, gamma) +
		(1.-eta) * transmission(EF + (eta-1.)*V, eps, gamma);
}
