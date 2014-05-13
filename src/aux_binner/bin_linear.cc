/**
 * \file bin_linear.cc
 * \brief Implements linear binning.
 *
 * \author Matthew G.\ Reuter
 * \date May 2014
 */

#include "bin_linear.h"

double BinLinear::gmask(const double g) {
	return g;
}

double BinLinear::invgmask(const double u) {
	return u;
}

double BinLinear::dudg(const double g) {
	return 1.;
}
