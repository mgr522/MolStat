/* This file is a part of MolStat, which is distributed under the Creative
   Commons Attribution-NonCommercial 4.0 International Public License.
   MolStat (c) 2014, Northwestern University. */
/**
 * \file sym_one_site_channel.cc
 * \brief Tight-binding channel with one site that couples symmetrically to
 *    both electrodes.
 *
 * \author Matthew G.\ Reuter
 * \date October 2014
 */

#include "sym_one_site_channel.h"
#include <cmath>

namespace molstat {
namespace transport {

const std::size_t SymOneSiteChannel::Index_EF = TransportJunction::Index_EF;
const std::size_t SymOneSiteChannel::Index_V = TransportJunction::Index_V;
const std::size_t SymOneSiteChannel::Index_epsilon = 2;
const std::size_t SymOneSiteChannel::Index_gamma = 3;
const std::size_t SymOneSiteChannel::Index_a = 4;

std::vector<std::string> SymOneSiteChannel::get_names() const {
	std::vector<std::string> ret(3);

	// subtract two because we expect 2 parameters from the TransportJunction
	ret[Index_epsilon - 2] = "epsilon";
	ret[Index_gamma - 2] = "gamma";
	ret[Index_a - 2] = "a";

	return ret;
}

double SymOneSiteChannel::transmission(const double e, const double V,
	const double eps, const double gamma, const double a) {

	return gamma*gamma / ((e - eps - a*V)*(e - eps - a*V) + gamma*gamma);
}

double SymOneSiteChannel::DiffG(const std::valarray<double> &params) const {
	// unpack the parameters
	const double &ef = params[Index_EF];
	const double &V = params[Index_V];
	const double &eps = params[Index_epsilon];
	const double &gamma = params[Index_gamma];
	const double &a = params[Index_a];
	
	return (0.5 - a) * transmission(ef+0.5*V, V, eps, gamma, a) +
		(0.5 + a) * transmission(ef-0.5*V, V, eps, gamma, a);
}

double SymOneSiteChannel::StaticG(const std::valarray<double> &params) const {
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