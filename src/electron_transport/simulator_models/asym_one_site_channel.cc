/* This file is a part of MolStat, which is distributed under the Creative
   Commons Attribution-NonCommercial 4.0 International Public License.

   (c) 2014 Northwestern University. */

/**
 * \file asym_one_site_channel.cc
 * \brief Tight-binding channel with one site that couples asymmetrically to
 *    both electrodes.
 *
 * \author Matthew G.\ Reuter
 * \date October 2014
 */

#include "asym_one_site_channel.h"
#include <cmath>

namespace molstat {
namespace transport {

const std::size_t AsymOneSiteChannel::Index_EF = TransportJunction::Index_EF;
const std::size_t AsymOneSiteChannel::Index_V = TransportJunction::Index_V;
const std::size_t AsymOneSiteChannel::Index_epsilon = 2;
const std::size_t AsymOneSiteChannel::Index_gammaL = 3;
const std::size_t AsymOneSiteChannel::Index_gammaR = 4;
const std::size_t AsymOneSiteChannel::Index_a = 5;

std::vector<std::string> AsymOneSiteChannel::get_names() const
{
	std::vector<std::string> ret(4);

	// subtract two because we expect 2 parameters from the TransportJunction
	ret[Index_epsilon - 2] = "epsilon";
	ret[Index_gammaL - 2] = "gammal";
	ret[Index_gammaR - 2] = "gammar";
	ret[Index_a - 2] = "a";

	return ret;
}

double AsymOneSiteChannel::transmission(const double e, const double V,
	const double eps, const double gammal, const double gammar,
	const double a)
{
	return 4.*gammal*gammar / (4.*(e - eps - a*V)*(e - eps - a*V) +
		(gammal + gammar)*(gammal + gammar));
}

double AsymOneSiteChannel::ECurrent(const std::valarray<double> &params) const
{
	// unpack the model parameters
	const double &ef = params[Index_EF];
	const double &V = params[Index_V];
	const double &eps = params[Index_epsilon];
	const double &gammal = params[Index_gammaL];
	const double &gammar = params[Index_gammaR];
	const double &a = params[Index_a];

	return 2. * TransportJunction::qc * gammal*gammar / (gammal + gammar) *
		(atan(2. * (ef-eps+(0.5-a)*V) / (gammal + gammar))
		- atan(2. * (ef-eps-(0.5+a)*V) / (gammal + gammar)));
}

double AsymOneSiteChannel::StaticG(const std::valarray<double> &params) const
{
	// unpack the model parameters
	const double &V = params[Index_V];

	return ECurrent(params) / (TransportJunction::qc * V);
}

double AsymOneSiteChannel::ZeroBiasG(const std::valarray<double> &params) const
{
	// unpack the model parameters
	const double &ef = params[Index_EF];
	const double &eps = params[Index_epsilon];
	const double &gammal = params[Index_gammaL];
	const double &gammar = params[Index_gammaR];
	const double &a = params[Index_a];

	return transmission(ef, 0., eps, gammal, gammar, a);
}

double AsymOneSiteChannel::DiffG(const std::valarray<double> &params) const
{
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
