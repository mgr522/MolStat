/* This file is a part of MolStat, which is distributed under the Creative
   Commons Attribution-NonCommercial 4.0 International Public License.
   MolStat (c) 2014, Northwestern University. */
/**
 * \file sym_interference.cc
 * \brief Tight-binding channel with two sites that have an interference
 *    feature.
 *
 * \author Matthew G.\ Reuter
 * \date February 2015
 */

#include "sym_interference.h"

namespace molstat {
namespace transport {

const std::size_t SymInterferenceChannel::Index_EF = TransportJunction::Index_EF;
const std::size_t SymInterferenceChannel::Index_V = TransportJunction::Index_V;
const std::size_t SymInterferenceChannel::Index_epsilon = 2;
const std::size_t SymInterferenceChannel::Index_gamma = 3;
const std::size_t SymInterferenceChannel::Index_beta = 4;

std::vector<std::string> SymInterferenceChannel::get_names() const
{
	std::vector<std::string> ret(3);

	// subtract two because we expect 2 parameters from the TransportJunction
	ret[Index_epsilon - 2] = "epsilon";
	ret[Index_gamma - 2] = "gamma";
	ret[Index_beta - 2] = "beta";

	return ret;
}

double SymInterferenceChannel::transmission(const double e, const double eps,
	const double gamma, const double beta)
{
	const double temp1 = e - eps;
	const double temp2 = temp1*temp1 - beta*beta;
	return gamma*gamma*temp1*temp1 / (temp2*temp2 + temp1*temp1*gamma*gamma);
}

double SymInterferenceChannel::ZeroBiasG(const std::valarray<double> &params) const
{
	// unpack the parameters
	const double &ef = params[Index_EF];
	const double &eps = params[Index_epsilon];
	const double &gamma = params[Index_gamma];
	const double &beta = params[Index_beta];
	
	return transmission(ef, eps, gamma, beta);
}

} // namespace molstat::transport
} // namespace molstat