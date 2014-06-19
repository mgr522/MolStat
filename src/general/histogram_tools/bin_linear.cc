/* This file is a part of MolStat, which is distributed under the Creative
   Commons Attribution-NonCommercial 4.0 International Public License.
   MolStat (c) 2014, Northwestern University. */
/**
 * \file bin_linear.cc
 * \brief Implements linear binning.
 *
 * \author Matthew G.\ Reuter
 * \date June 2014
 */

#include "bin_linear.h"

double BinLinear::mask(const double x) const {
	return x;
}

double BinLinear::invmask(const double u) const {
	return u;
}

double BinLinear::dmaskdx(const double x) const {
	return 1.;
}
