/**
 * \file aux_simulator/symmetric_voltage_two_site.cc
 * \brief Implementation of the symmetric-coupling, voltage-dependent two-site
 *    tight-binding model for calculating conductances.
 *
 * \author Matthew G.\ Reuter
 * \date May 2014
 */

#include "symmetric_voltage_two_site.h"
#include <cmath>
#include <complex>
#include <string>
#include <vector>
#include <stdexcept>

using namespace std;

SymmetricVoltageTwoSiteModel::SymmetricVoltageTwoSiteModel(
	const shared_ptr<const RandomDistribution> &eps,
	const shared_ptr<const RandomDistribution> &gamma,
	const shared_ptr<const RandomDistribution> &beta)
	: dist_eps(eps), dist_gamma(gamma), dist_beta(beta) {}

double SymmetricVoltageTwoSiteModel::static_conductance(
	shared_ptr<gsl_rng> r, const double EF, const double eta,
	const double V) const {

	// get model parameters from the random distributions
	double eps = dist_eps->sample(r);
	double gamma = dist_gamma->sample(r);
	double beta = dist_beta->sample(r);

	return SymmetricVoltageTwoSiteModel::static_conductance(EF, V, eta, eps,
		gamma, beta);
}

double SymmetricVoltageTwoSiteModel::diff_conductance(
	shared_ptr<gsl_rng> r, const double EF, const double eta,
	const double V) const {

	// get model parameters from the random distributions
	double eps = dist_eps->sample(r);
	double gamma = dist_gamma->sample(r);
	double beta = dist_beta->sample(r);

	return SymmetricVoltageTwoSiteModel::diff_conductance(EF, V, eta, eps,
		gamma, beta);
}

double SymmetricVoltageTwoSiteModel::zero_bias_conductance(
	shared_ptr<gsl_rng> r, const double EF) const {

	// get model parameters from the random distributions
	double eps = dist_eps->sample(r);
	double gamma = dist_gamma->sample(r);
	double beta = dist_beta->sample(r);

	return SymmetricVoltageTwoSiteModel::transmission(EF, 0., eps, gamma, beta);
}

double SymmetricVoltageTwoSiteModel::transmission(const double E,
	const double V, const double eps, const double gamma, const double beta) {

	double temp = 4.*(E-eps)*(E-eps) - 4.*beta*beta - gamma*gamma - V*V;

	return 16.*gamma*gamma*beta*beta /
		(temp*temp + 16.*gamma*gamma*(E-eps)*(E-eps));
}

double SymmetricVoltageTwoSiteModel::static_c_integral(const double z,
	const double V, const double eps, const double gamma, const double beta) {

	complex<double> arg(gamma, -sqrt(4.*beta*beta + V*V));

	return 4.*beta*beta*gamma / sqrt(4.*beta*beta + V*V) *
		imag(atan(2.*(z-eps)/arg)/arg);

}

double SymmetricVoltageTwoSiteModel::static_conductance(const double EF,
	const double V, const double eta, const double eps, const double gamma,
	const double beta) {

	return (SymmetricVoltageTwoSiteModel::static_c_integral(EF+eta*V, V,
		eps, gamma, beta) - SymmetricVoltageTwoSiteModel::static_c_integral(
		EF+(eta-1.)*V, V, eps, gamma, beta)) / V;
}

double SymmetricVoltageTwoSiteModel::diff_c_integral(const double z,
	const double V, const double gamma, const double beta) {

	const double bv = 4.*beta*beta + V*V;
	const double bvg = bv + gamma*gamma;
	const complex<double> arctan = atan(2.*z /
		complex<double>(gamma, -sqrt(bv)));

	return 8.*V*gamma*gamma*beta*beta*z*(4.*z*z + gamma*gamma - 3.*bv) /
		(bv*bvg*(16.*z*z*z*z + 8.*(gamma*gamma - bv)*z*z + bvg*bvg))

		- 8.*V*gamma*beta*beta / (bvg*bvg) * real(arctan)

		- 4.*V*gamma*gamma*beta*beta*(3.*bv + gamma*gamma) /
			(bvg*bvg*bv*sqrt(bv)) * imag(arctan);
}

double SymmetricVoltageTwoSiteModel::diff_conductance(const double EF,
	const double V, const double eta, const double eps, const double gamma,
	const double beta) {

	return eta * SymmetricVoltageTwoSiteModel::transmission(EF + eta*V, V,
		eps, gamma, beta)

		+ (1.-eta) * SymmetricVoltageTwoSiteModel::transmission(EF + (eta-1.)*V,
		V, eps, gamma, beta)

		+ SymmetricVoltageTwoSiteModel::diff_c_integral(EF - eps + eta*V, V,
		gamma, beta)

		- SymmetricVoltageTwoSiteModel::diff_c_integral(EF - eps + (eta-1.)*V, V,
		gamma, beta);
}
