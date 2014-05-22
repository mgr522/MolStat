/* This file is a part of MolStat, which is distributed under the Creative
   Commons Attribution-NonCommercial 4.0 International Public License.
   MolStat (c) 2014, Northwestern University. */
/**
 * \file constant.cc
 * \brief Implementation of the constant distribution.
 *
 * \author Matthew G.\ Reuter
 * \date May 2014
 */

#include <gsl/gsl_randist.h>
#include "constant.h"

ConstantDistribution::ConstantDistribution(const double val)
	: RandomDistribution(), value(val) {}

double ConstantDistribution::sample(shared_ptr<gsl_rng> r) const {
	return value;
}