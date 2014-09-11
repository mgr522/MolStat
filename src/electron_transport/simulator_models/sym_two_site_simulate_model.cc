/* This file is a part of MolStat, which is distributed under the Creative
   Commons Attribution-NonCommercial 4.0 International Public License.
   MolStat (c) 2014, Northwestern University. */
/**
 * \file sym_two_site_simulate_model.cc
 * \brief Tight-binding model of a two-site chain that couples symmetrically to
 *    both electrodes.
 *
 * \author Matthew G.\ Reuter
 * \date September 2014
 */

#include "sym_two_site_simulate_model.h"
#include <cmath>
#include <complex>

using namespace std;

// if the order of the following list is changed, the unpack_parameters
// function MUST also be updated
const vector<string> SymTwoSiteSimulateModel::parameters =
	{"ef", "v", "epsilon", "gamma", "beta"};

void SymTwoSiteSimulateModel::unpack_parameters(const std::vector<double> &vec,
	double &ef, double &v, double &epsilon, double &gamma, double &beta) {

	ef = vec[0];
	v = vec[1];
	epsilon = vec[2];
	gamma = vec[3];
	beta = vec[4];
}

SymTwoSiteSimulateModel::SymTwoSiteSimulateModel(
	const std::map<std::string, shared_ptr<RandomDistribution>> &avail)
	: SimulateModel(avail, parameters) {
}

std::array<double, 2> SymTwoSiteSimulateModel::DiffG(shared_ptr<gsl_rng> r)
	const {

	vector<double> params(5);
	double ef, v, eps, gamma, beta;

	// generate and unpack the parameters
	sample(r, params);
	unpack_parameters(params, ef, v, eps, gamma, beta);

	return {{ v, diff_conductance(params) }};
}

std::array<double, 2> SymTwoSiteSimulateModel::StaticG(shared_ptr<gsl_rng> r)
	const {

	vector<double> params(5);
	double ef, v, eps, gamma, beta;

	// generate and unpack the parameters
	sample(r, params);
	unpack_parameters(params, ef, v, eps, gamma, beta);

	return {{ v, static_conductance(params) }};
}

double SymTwoSiteSimulateModel::transmission(const double e, const double v,
	const double eps, const double gamma, const double beta) {

	double temp = 4.*(e-eps)*(e-eps) - 4.*beta*beta - gamma*gamma;

	return 16.*gamma*gamma*beta*beta /
		(temp*temp + 16.*gamma*gamma*(e-eps)*(e-eps));
}

double SymTwoSiteSimulateModel::static_c_integral(const double z,
	const double eps, const double gamma, const double beta) {

	return 2.*beta*gamma / (4.*beta*beta + gamma*gamma) *
		real(complex<double>(gamma, 2.*beta)
		* atanh(2.*(z-eps) / complex<double>(2.*beta, gamma)));
}

double SymTwoSiteSimulateModel::static_conductance(
	const std::vector<double> &vec) {

	double ef, v, eps, gamma, beta;

	// unpack the model parameters
	unpack_parameters(vec, ef, v, eps, gamma, beta);

	return (static_c_integral(ef + 0.5*v, eps, gamma, beta) -
		static_c_integral(ef - 0.5*v, eps, gamma, beta)) / v;
}

double SymTwoSiteSimulateModel::diff_conductance(
	const std::vector<double> &vec) {

	double ef, v, eps, gamma, beta;

	// unpack the parameters
	unpack_parameters(vec, ef, v, eps, gamma, beta);

	return 0.5*transmission(ef + 0.5*v, v, eps, gamma, beta) +
		0.5*transmission(ef - 0.5*v, v, eps, gamma, beta);
}
