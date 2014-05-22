/* This file is a part of MolStat, which is distributed under the Creative
   Commons Attribution-NonCommercial 4.0 International Public License.
   MolStat (c) 2014, Northwestern University. */
/**
 * \file simulator_models/symmetric_one_site.cc
 * \brief Implementation of the symmetric-coupling, one-site model for
 *        calculating conductances.
 *
 * \author Matthew G.\ Reuter
 * \date April 2014
 */

#include "symmetric_one_site.h"
#include <cmath>
#include <string>
#include <vector>
#include <stdexcept>

using namespace std;

SymmetricOneSiteModel::SymmetricOneSiteModel(
	const shared_ptr<const RandomDistribution> &eps,
	const shared_ptr<const RandomDistribution> &gamma)
	: dist_eps(eps), dist_gamma(gamma) {}

double SymmetricOneSiteModel::static_conductance(
	shared_ptr<gsl_rng> r, const double EF, const double eta,
	const double V) const {

	// get model parameters from the random distributions
	double eps = dist_eps->sample(r);
	double gamma = dist_gamma->sample(r);

	return static_conductance(EF, V, eta, eps, gamma);
}

double SymmetricOneSiteModel::diff_conductance(
	shared_ptr<gsl_rng> r, const double EF, const double eta,
	const double V) const {

	// get model parameters from the random distributions
	double eps = dist_eps->sample(r);
	double gamma = dist_gamma->sample(r);

	return diff_conductance(EF, V, eta, eps, gamma);
}

double SymmetricOneSiteModel::zero_bias_conductance(
	shared_ptr<gsl_rng> r, const double EF) const {

	// get model parameters from the random distributions
	double eps = dist_eps->sample(r);
	double gamma = dist_gamma->sample(r);

	return transmission(EF, eps, gamma);
}

double SymmetricOneSiteModel::transmission(const double E,
	const double eps, const double gamma) {

	return gamma*gamma / ((E-eps)*(E-eps) + gamma*gamma);
}

double SymmetricOneSiteModel::static_conductance(const double EF,
	const double V, const double eta, const double eps, const double gamma) {

	return gamma / V
		* (atan((EF-eps+eta*V) / gamma) - atan((EF-eps+(eta-1.)*V) / gamma));
}

double SymmetricOneSiteModel::diff_conductance(const double EF,
	const double V, const double eta, const double eps, const double gamma) {

	return eta * transmission(EF + eta*V, eps, gamma) +
		(1.-eta) * transmission(EF + (eta-1.)*V, eps, gamma);
}
