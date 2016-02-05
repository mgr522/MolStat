/* This file is a part of MolStat, which is distributed under the Creative
   Commons Attribution-NonCommercial 4.0 International Public License.

   (c) 2014 Northwestern University. */

/**
 * \file asym_two_site_channel.cc
 * \brief Tight-binding channel with a two-site chain that couples
 *    asymmetrically to both electrodes.
 *
 * \author Matthew G.\ Reuter
 * \date October 2014
 */

#include "asym_two_site_channel.h"
#include <cmath>
#include <complex>

namespace molstat {
namespace transport {

const std::size_t AsymTwoSiteChannel::Index_EF = TransportJunction::Index_EF;
const std::size_t AsymTwoSiteChannel::Index_V = TransportJunction::Index_V;
const std::size_t AsymTwoSiteChannel::Index_epsilon = 2;
const std::size_t AsymTwoSiteChannel::Index_gammaL = 3;
const std::size_t AsymTwoSiteChannel::Index_gammaR = 4;
const std::size_t AsymTwoSiteChannel::Index_beta = 5;

std::vector<std::string> AsymTwoSiteChannel::get_names() const
{
	std::vector<std::string> ret(4);

	// subtract two because we expect 2 parameters from the TranportJunction
	ret[Index_epsilon - 2] = "epsilon";
	ret[Index_gammaL - 2] = "gammal";
	ret[Index_gammaR - 2] = "gammar";
	ret[Index_beta - 2] = "beta";

	return ret;
}

double AsymTwoSiteChannel::transmission(const double e, const double V,
	const double eps, const double gammal, const double gammar,
	const double beta)
{
	double temp = 4.*(e-eps)*(e-eps) - 4.*beta*beta - gammal*gammar;

	return 16.*gammal*gammar*beta*beta /
		(temp*temp + 4.*(gammal+gammar)*(gammal+gammar)*(e-eps)*(e-eps));
}

double AsymTwoSiteChannel::current_integral(const double z,
	const double eps, const double gammal, const double gammar,
	const double beta)
{
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

double AsymTwoSiteChannel::ECurrent(const std::valarray<double> &params) const
{
	// unpack the parameters
	const double &ef = params[Index_EF];
	const double &V = params[Index_V];
	const double &eps = params[Index_epsilon];
	const double &gammal = params[Index_gammaL];
	const double &gammar = params[Index_gammaR];
	const double &beta = params[Index_beta];

	return TransportJunction::qc *
		(current_integral(ef + 0.5*V, eps, gammal, gammar, beta) -
		 current_integral(ef - 0.5*V, eps, gammal, gammar, beta));
}

double AsymTwoSiteChannel::StaticG(const std::valarray<double> &params) const
{
	// unpack the parameters
	const double &V = params[Index_V];

	return ECurrent(params) / (TransportJunction::qc * V);
}

double AsymTwoSiteChannel::ZeroBiasG(const std::valarray<double> &params) const
{
	// unpack the parameters
	const double &ef = params[Index_EF];
	const double &eps = params[Index_epsilon];
	const double &gammal = params[Index_gammaL];
	const double &gammar = params[Index_gammaR];
	const double &beta = params[Index_beta];

	return transmission(ef, 0., eps, gammal, gammar, beta);
}

double AsymTwoSiteChannel::DiffG(const std::valarray<double> &params) const
{
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
