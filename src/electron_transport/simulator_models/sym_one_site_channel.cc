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
namespace transport {

const std::size_t SymOneSiteSimulateModel::Index_EF = 0;
const std::size_t SymOneSiteSimulateModel::Index_V = 1;
const std::size_t SymOneSiteSimulateModel::Index_epsilon = 2;
const std::size_t SymOneSiteSimulateModel::Index_gamma = 3;
const std::size_t SymOneSiteSimulateModel::Index_a = 4;

SymOneSiteSimulateModel::SymOneSiteSimulateModel(
	const std::map<std::string, std::shared_ptr<RandomDistribution>> &avail)
	: SimulateModel(avail,
	                order_from_map({{Index_EF, "ef"},
	                                {Index_V, "v"},
	                                {Index_epsilon, "epsilon"},
	                                {Index_gamma, "gamma"},
	                                {Index_a, "a"}})) {
}

double SymOneSiteSimulateModel::transmission(const double e, const double V,
	const double eps, const double gamma, const double a) {

	return gamma*gamma / ((e - eps - a*V)*(e - eps - a*V) + gamma*gamma);
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

} // namespace molstat::transport
} // namespace molstat