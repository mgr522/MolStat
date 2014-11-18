/* This file is a part of MolStat, which is distributed under the Creative
   Commons Attribution-NonCommercial 4.0 International Public License.
   MolStat (c) 2014, Northwestern University. */
/**
 * \file sym_two_site_channel.cc
 * \brief Tight-binding channel of a two-site chain that couples symmetrically
 *    to both electrodes.
 *
 * \author Matthew G.\ Reuter
 * \date October 2014
 */

#include "sym_two_site_channel.h"
#include <cmath>
#include <complex>

namespace molstat {
namespace transport {

const std::size_t SymTwoSiteChannel::Index_EF = TransportJunction::Index_EF;
const std::size_t SymTwoSiteChannel::Index_V = TransportJunction::Index_V;
const std::size_t SymTwoSiteChannel::Index_epsilon = 2;
const std::size_t SymTwoSiteChannel::Index_gamma = 3;
const std::size_t SymTwoSiteChannel::Index_beta = 4;

std::vector<std::string> SymTwoSiteChannel::get_names() const
{
	std::vector<std::string> ret(3);

	// subtract two because we expect 2 parameters from the TransportJunction
	ret[Index_epsilon - 2] = "epsilon";
	ret[Index_gamma - 2] = "gamma";
	ret[Index_beta - 2] = "beta";

	return ret;
}

double SymTwoSiteChannel::transmission(const double e, const double V,
	const double eps, const double gamma, const double beta)
{
	double temp = 4.*(e-eps)*(e-eps) - 4.*beta*beta - gamma*gamma;

	return 16.*gamma*gamma*beta*beta /
		(temp*temp + 16.*gamma*gamma*(e-eps)*(e-eps));
}

double SymTwoSiteChannel::static_c_integral(const double z,
	const double eps, const double gamma, const double beta)
{
	return 2.*beta*gamma / (4.*beta*beta + gamma*gamma) *
		real(std::complex<double>(gamma, 2.*beta)
		* atanh(2.*(z-eps) / std::complex<double>(2.*beta, gamma)));
}

double SymTwoSiteChannel::StaticG(const std::valarray<double> &params) const
{
	// unpack the parameters
	const double &ef = params[Index_EF];
	const double &V = params[Index_V];
	const double &eps = params[Index_epsilon];
	const double &gamma = params[Index_gamma];
	const double &beta = params[Index_beta];
	
	return (static_c_integral(ef + 0.5*V, eps, gamma, beta) -
		static_c_integral(ef - 0.5*V, eps, gamma, beta)) / V;
}

double SymTwoSiteChannel::DiffG(const std::valarray<double> &params) const
{
	// unpack the parameters
	const double &ef = params[Index_EF];
	const double &V = params[Index_V];
	const double &eps = params[Index_epsilon];
	const double &gamma = params[Index_gamma];
	const double &beta = params[Index_beta];
	
	return 0.5*transmission(ef + 0.5*V, V, eps, gamma, beta) +
		0.5*transmission(ef - 0.5*V, V, eps, gamma, beta);
}

} // namespace molstat::transport
} // namespace molstat
