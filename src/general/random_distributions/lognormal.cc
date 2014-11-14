/* This file is a part of MolStat, which is distributed under the Creative
   Commons Attribution-NonCommercial 4.0 International Public License.
   MolStat (c) 2014, Northwestern University. */
/**
 * \internal
 * \file lognormal.cc
 * \brief Implementation of the lognormal distribution.
 *
 * \author Matthew G.\ Reuter
 * \date May 2014
 * \endinternal
 */

#include <gsl/gsl_randist.h>
#include "lognormal.h"

namespace molstat {

LognormalDistribution::LognormalDistribution(const double zeta_,
	const double sigma_)
	: RandomDistribution(), zeta(zeta_), sigma(sigma_)
{
}

double LognormalDistribution::sample(gsl_rng_ptr &r) const
{
	return gsl_ran_lognormal(r.get(), zeta, sigma);
}

std::string LognormalDistribution::info() const
{
	return "Lognormal: mean = " + std::to_string(zeta) +
		" and stdev = " + std::to_string(sigma) + " (log space).";
}

} // namespace molstat
