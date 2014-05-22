/* This file is a part of MolStat, which is distributed under the Creative
   Commons Attribution-NonCommercial 4.0 International Public License.
   MolStat (c) 2014, Northwestern University. */
/**
 * \file simulator_models/symmetric_two_site.cc
 * \brief Implementation of the symmetric-coupling, two-site tight-binding
 *    model for calculating conductances.
 *
 * \author Matthew G.\ Reuter
 * \date May 2014
 */

#include "symmetric_two_site.h"
#include <cmath>
#include <complex>
#include <string>
#include <vector>
#include <stdexcept>

using namespace std;

SymmetricTwoSiteModel::SymmetricTwoSiteModel(
	const shared_ptr<const RandomDistribution> &eps,
	const shared_ptr<const RandomDistribution> &gamma,
	const shared_ptr<const RandomDistribution> &beta)
	: dist_eps(eps), dist_gamma(gamma), dist_beta(beta) {}

double SymmetricTwoSiteModel::static_conductance(
	shared_ptr<gsl_rng> r, const double EF, const double eta,
	const double V) const {

	// get model parameters from the random distributions
	double eps = dist_eps->sample(r);
	double gamma = dist_gamma->sample(r);
	double beta = dist_beta->sample(r);

	return SymmetricTwoSiteModel::static_conductance(EF, V, eta, eps,
		gamma, beta);
}

double SymmetricTwoSiteModel::diff_conductance(
	shared_ptr<gsl_rng> r, const double EF, const double eta,
	const double V) const {

	// get model parameters from the random distributions
	double eps = dist_eps->sample(r);
	double gamma = dist_gamma->sample(r);
	double beta = dist_beta->sample(r);

	return SymmetricTwoSiteModel::diff_conductance(EF, V, eta, eps,
		gamma, beta);
}

double SymmetricTwoSiteModel::zero_bias_conductance(
	shared_ptr<gsl_rng> r, const double EF) const {

	// get model parameters from the random distributions
	double eps = dist_eps->sample(r);
	double gamma = dist_gamma->sample(r);
	double beta = dist_beta->sample(r);

	return SymmetricTwoSiteModel::transmission(EF, eps, gamma, beta);
}

double SymmetricTwoSiteModel::transmission(const double E,
	const double eps, const double gamma, const double beta) {

	double temp = 4.*(E-eps)*(E-eps) - 4.*beta*beta - gamma*gamma;

	return 16.*gamma*gamma*beta*beta /
		(temp*temp + 16.*gamma*gamma*(E-eps)*(E-eps));
}

double SymmetricTwoSiteModel::static_c_integral(const double z,
	const double eps, const double gamma, const double beta) {

	return 2.*beta*gamma / (4.*beta*beta + gamma*gamma) *
		real(complex<double>(gamma, 2.*beta)
		* atanh(2.*(z-eps) / complex<double>(2.*beta, gamma)));

}

double SymmetricTwoSiteModel::static_conductance(const double EF,
	const double V, const double eta, const double eps, const double gamma,
	const double beta) {

	return (SymmetricTwoSiteModel::static_c_integral(EF+eta*V, eps, gamma,
		beta) - SymmetricTwoSiteModel::static_c_integral(EF+(eta-1.)*V, eps,
		gamma, beta)) / V;
}

double SymmetricTwoSiteModel::diff_conductance(const double EF,
	const double V, const double eta, const double eps, const double gamma,
	const double beta) {

	return eta * SymmetricTwoSiteModel::transmission(EF + eta*V, eps, gamma,
		beta)

		+ (1.-eta) * SymmetricTwoSiteModel::transmission(EF + (eta-1.)*V, eps,
		gamma, beta);
}
