/* This file is a part of MolStat, which is distributed under the Creative
   Commons Attribution-NonCommercial 4.0 International Public License.
   MolStat (c) 2014, Northwestern University. */
/**
 * \file bin_linear.cc
 * \brief Implements linear binning.
 *
 * \author Matthew G.\ Reuter
 * \date May 2014
 */

#include "bin_linear.h"

double BinLinear::gmask(const double g) const {
	return g;
}

double BinLinear::invgmask(const double u) const {
	return u;
}

double BinLinear::dudg(const double g) const {
	return 1.;
}
