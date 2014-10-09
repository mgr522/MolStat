/* This file is a part of MolStat, which is distributed under the Creative
   Commons Attribution-NonCommercial 4.0 International Public License.
   MolStat (c) 2014, Northwestern University. */
/**
 * \file sym_one_site_simulate_model.cc
 * \brief Tight-binding model with one site that couples symmetrically to
 *    both electrodes.
 *
 * \author Matthew G.\ Reuter
 * \date October 2014
 */

#include "sym_one_site_simulate_model.h"
#include <cmath>

namespace molstat {

// set up the ordering of parameters
const std::map<std::size_t, std::string> SymOneSiteSimulateModel::param_order =
	{{Index_EF, "ef"}, {Index_V, "v"}, {Index_epsilon, "epsilon"},
	 {Index_gamma, "gamma"}, {Index_a, "a"}};

SymOneSiteSimulateModel::SymOneSiteSimulateModel(
	const std::map<std::string, std::shared_ptr<RandomDistribution>> &avail)
	: SimulateModel(avail, order_from_map(param_order)) {
}

double SymOneSiteSimulateModel::transmission(const double e, const double v,
	const double eps, const double gamma, const double a) {

	return gamma*gamma / ((e - eps - a*v)*(e - eps - a*v) + gamma*gamma);
}

double SymOneSiteSimulateModel::AppBias(const std::array<double, 5> &params)
	const {

	return params[Index_V];
}

double SymOneSiteSimulateModel::DiffG(const std::array<double, 5> &params)
	const {

	// unpack the parameters
	const double &ef = params[Index_EF];
	const double &V = params[Index_V];
	const double &eps = params[Index_epsilon];
	const double &gamma = params[Index_gamma];
	const double &a = params[Index_a];
	
	return (0.5 - a) * transmission(ef+0.5*V, V, eps, gamma, a) +
		(0.5 + a) * transmission(ef-0.5*V, V, eps, gamma, a);
}

double SymOneSiteSimulateModel::StaticG(const std::array<double, 5> &params)
	const {

	// unpack the parameters
	const double &ef = params[Index_EF];
	const double &V = params[Index_V];
	const double &eps = params[Index_epsilon];
	const double &gamma = params[Index_gamma];
	const double &a = params[Index_a];

	return gamma / V *
		(atan((ef-eps+(0.5-a)*V) / gamma) - atan((ef-eps-(0.5+a)*V) / gamma));
}

} // namespace molstat