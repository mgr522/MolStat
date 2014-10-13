/* This file is a part of MolStat, which is distributed under the Creative
   Commons Attribution-NonCommercial 4.0 International Public License.
   MolStat (c) 2014, Northwestern University. */
/**
 * \file asym_two_site_simulate_model.cc
 * \brief Tight-binding model of a two-site chain that couples asymmetrically
 *    to both electrodes.
 *
 * \author Matthew G.\ Reuter
 * \date October 2014
 */

#include "asym_two_site_simulate_model.h"
#include <cmath>
#include <complex>

namespace molstat {
namespace transport {

const std::size_t AsymTwoSiteSimulateModel::Index_EF = 0;
const std::size_t AsymTwoSiteSimulateModel::Index_V = 1;
const std::size_t AsymTwoSiteSimulateModel::Index_epsilon = 2;
const std::size_t AsymTwoSiteSimulateModel::Index_gammaL = 3;
const std::size_t AsymTwoSiteSimulateModel::Index_gammaR = 4;
const std::size_t AsymTwoSiteSimulateModel::Index_beta = 5;

AsymTwoSiteSimulateModel::AsymTwoSiteSimulateModel(
	const std::map<std::string, std::shared_ptr<RandomDistribution>> &avail)
	: SimulateModel(avail,
	                order_from_map({{Index_EF, "ef"},
	                                {Index_V, "v"},
	                                {Index_epsilon, "epsilon"},
	                                {Index_gammaL, "gammal"},
	                                {Index_gammaR, "gammar"},
	                                {Index_beta, "beta"}})) {
}

double AsymTwoSiteSimulateModel::transmission(const double e, const double V,
	const double eps, const double gammal, const double gammar,
	const double beta) {

	double temp = 4.*(e-eps)*(e-eps) - 4.*beta*beta - gammal*gammar;

	return 16.*gammal*gammar*beta*beta /
		(temp*temp + 4.*(gammal+gammar)*(gammal+gammar)*(e-eps)*(e-eps));
}

double AsymTwoSiteSimulateModel::AppBias(const std::array<double, 6> &params)
	const {

	return params[Index_V];
}

double AsymTwoSiteSimulateModel::static_c_integral(const double z,
	const double eps, const double gammal, const double gammar,
	const double beta) {

	std::complex<double> bgg = sqrt(std::complex<double>((gammal-gammar)*(gammal-gammar)
		- 16.*beta*beta, 0.));

	std::complex<double> denom1 = sqrt(-8.*beta*beta + gammal*gammal + gammar*gammar
		- (gammal + gammar) * bgg);

	std::complex<double> denom2 = sqrt(-8.*beta*beta + gammal*gammal + gammar*gammar
		+ (gammal + gammar) * bgg);

	return sqrt(128.)*gammal*gammar*beta*beta / (gammal + gammar) * real(
		(atan(sqrt(8.)*(z-eps) / denom1) / denom1 -
		 atan(sqrt(8.)*(z-eps) / denom2) / denom2) / bgg
		);
}

double AsymTwoSiteSimulateModel::StaticG(const std::array<double, 6> &params)
	const {

	// unpack the parameters
	const double &ef = params[Index_EF];
	const double &V = params[Index_V];
	const double &eps = params[Index_epsilon];
	const double &gammal = params[Index_gammaL];
	const double &gammar = params[Index_gammaR];
	const double &beta = params[Index_beta];

	return (static_c_integral(ef + 0.5*V, eps, gammal, gammar, beta) -
		static_c_integral(ef - 0.5*V, eps, gammal, gammar, beta)) / V;
}

double AsymTwoSiteSimulateModel::DiffG(const std::array<double, 6> &params)
	const {

	// unpack the parameters
	const double &ef = params[Index_EF];
	const double &V = params[Index_V];
	const double &eps = params[Index_epsilon];
	const double &gammal = params[Index_gammaL];
	const double &gammar = params[Index_gammaR];
	const double &beta = params[Index_beta];

	return 0.5*transmission(ef + 0.5*V, V, eps, gammal, gammar, beta) +
		0.5*transmission(ef - 0.5*V, V, eps, gammal, gammar, beta);
}

} // namespace molstat::transport
} // namespace molstat