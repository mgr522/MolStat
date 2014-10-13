/* This file is a part of MolStat, which is distributed under the Creative
   Commons Attribution-NonCommercial 4.0 International Public License.
   MolStat (c) 2014, Northwestern University. */
/**
 * \file asym_one_site_simulate_model.cc
 * \brief Tight-binding model with one site that couples asymmetrically to
 *    both electrodes.
 *
 * \author Matthew G.\ Reuter
 * \date October 2014
 */

#include "asym_one_site_simulate_model.h"
#include <cmath>

namespace molstat {
namespace transport {

const std::size_t AsymOneSiteSimulateModel::Index_EF = 0;
const std::size_t AsymOneSiteSimulateModel::Index_V = 1;
const std::size_t AsymOneSiteSimulateModel::Index_epsilon = 2;
const std::size_t AsymOneSiteSimulateModel::Index_gammaL = 3;
const std::size_t AsymOneSiteSimulateModel::Index_gammaR = 4;
const std::size_t AsymOneSiteSimulateModel::Index_a = 5;

AsymOneSiteSimulateModel::AsymOneSiteSimulateModel(
	const std::map<std::string, std::shared_ptr<RandomDistribution>> &avail)
	: SimulateModel(avail,
	                order_from_map({{Index_EF, "ef"},
	                                {Index_V, "v"},
	                                {Index_epsilon, "epsilon"},
	                                {Index_gammaL, "gammal"},
	                                {Index_gammaR, "gammar"},
	                                {Index_a, "a"}})) {
}

double AsymOneSiteSimulateModel::transmission(const double e, const double V,
	const double eps, const double gammal, const double gammar,
	const double a) {

	return 4.*gammal*gammar / (4.*(e - eps - a*V)*(e - eps - a*V) +
		(gammal + gammar)*(gammal + gammar));
}

double AsymOneSiteSimulateModel::AppBias(const std::array<double, 6> &params)
	const {

	return params[Index_V];
}

double AsymOneSiteSimulateModel::StaticG(const std::array<double, 6> &params
	) const {

	// unpack the model parameters
	const double &ef = params[Index_EF];
	const double &V = params[Index_V];
	const double &eps = params[Index_epsilon];
	const double &gammal = params[Index_gammaL];
	const double &gammar = params[Index_gammaR];
	const double &a = params[Index_a];

	return 2.*gammal*gammar / (V*(gammal + gammar)) *
		(atan(2. * (ef-eps+(0.5-a)*V) / (gammal + gammar))
		- atan(2. * (ef-eps-(0.5+a)*V) / (gammal + gammar)));
}

double AsymOneSiteSimulateModel::DiffG(const std::array<double, 6> &params
	) const {

	// unpack the model parameters
	const double &ef = params[Index_EF];
	const double &V = params[Index_V];
	const double &eps = params[Index_epsilon];
	const double &gammal = params[Index_gammaL];
	const double &gammar = params[Index_gammaR];
	const double &a = params[Index_a];

	return (0.5 - a) * transmission(ef+0.5*V, V, eps, gammal, gammar, a) +
		(0.5 + a) * transmission(ef-0.5*V, V, eps, gammal, gammar, a);
}

} // namespace molstat::transport
} // namespace molstat
