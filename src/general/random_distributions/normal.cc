/* This file is a part of MolStat, which is distributed under the Creative
   Commons Attribution-NonCommercial 4.0 International Public License.
   MolStat (c) 2014, Northwestern University. */
/**
 * \internal
 * \file normal.cc
 * \brief Implementation of the normal distribution.
 *
 * \author Matthew G.\ Reuter
 * \date May 2014
 * \endinternal
 */

#include <gsl/gsl_randist.h>
#include "normal.h"

NormalDistribution::NormalDistribution(const double mean_, const double stdev_)
	: RandomDistribution(), mean(mean_), stdev(stdev_) {}

double NormalDistribution::sample(gsl_rng_ptr &r) const {
	return gsl_ran_gaussian(r.get(), stdev) + mean;
}

std::string NormalDistribution::info() const {
	return "Normal distribution with mean " + std::to_string(mean) +
		" and standard deviation " + std::to_string(stdev) + ".";
}
