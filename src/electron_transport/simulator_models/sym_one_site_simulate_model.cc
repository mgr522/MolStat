/* This file is a part of MolStat, which is distributed under the Creative
   Commons Attribution-NonCommercial 4.0 International Public License.
   MolStat (c) 2014, Northwestern University. */
/**
 * \file sym_one_site_simulate_model.cc
 * \brief Tight-binding model with one site that couples symmetrically to
 *    both electrodes.
 *
 * \author Matthew G.\ Reuter
 * \date September 2014
 */

#include "sym_one_site_simulate_model.h"
#include <cmath>

namespace molstat {

// if the order of the following list is changed, the unpack_parameters
// function MUST also be updated
const std::vector<std::string> SymOneSiteSimulateModel::parameters =
	{"ef", "v", "epsilon", "gamma", "a"};

void SymOneSiteSimulateModel::unpack_parameters(const std::vector<double> &vec,
	double &ef, double &v, double &epsilon, double &gamma, double &a) {

	ef = vec[0];
	v = vec[1];
	epsilon = vec[2];
	gamma = vec[3];
	a = vec[4];
}

SymOneSiteSimulateModel::SymOneSiteSimulateModel(
	const std::map<std::string, std::shared_ptr<RandomDistribution>> &avail)
	: SimulateModel(avail, parameters) {
}

std::array<double, 2> SymOneSiteSimulateModel::DiffG(gsl_rng_ptr &r) const {
	std::vector<double> params(5);
	double ef, v, eps, gamma, a;

	// generate and unpack the parameters
	sample(r, params);
	unpack_parameters(params, ef, v, eps, gamma, a);

	return {{ v, diff_conductance(params) }};
}

std::array<double, 2> SymOneSiteSimulateModel::StaticG(gsl_rng_ptr &r) const {
	std::vector<double> params(5);
	double ef, v, eps, gamma, a;

	// generate and unpack the parameters
	sample(r, params);
	unpack_parameters(params, ef, v, eps, gamma, a);

	return {{ v, static_conductance(params) }};
}

double SymOneSiteSimulateModel::transmission(const double e, const double v,
	const double eps, const double gamma, const double a) {

	return gamma*gamma / ((e - eps - a*v)*(e - eps - a*v) + gamma*gamma);
}

double SymOneSiteSimulateModel::static_conductance(
	const std::vector<double> &vec) {

	double ef, v, eps, gamma, a;

	// unpack the model parameters
	unpack_parameters(vec, ef, v, eps, gamma, a);

	return gamma / v *
		(atan((ef-eps+(0.5-a)*v) / gamma) - atan((ef-eps-(0.5+a)*v) / gamma));
}

double SymOneSiteSimulateModel::diff_conductance(
	const std::vector<double> &vec) {

	double ef, v, eps, gamma, a;

	// unpack the parameters
	unpack_parameters(vec, ef, v, eps, gamma, a);

	return (0.5 - a) * transmission(ef+0.5*v, v, eps, gamma, a) +
		(0.5 + a) * transmission(ef-0.5*v, v, eps, gamma, a);
}

} // namespace molstat