/**
 * \file bin_log.cc
 * \brief Implements logarithmic binning.
 *
 * \author Matthew G.\ Reuter
 * \date May 2014
 */

#include "bin_log.h"
#include <cmath>

double BinLog::gmask(const double g) const {
	return log10(g);
}

double BinLog::invgmask(const double u) const {
	return pow(10., u);
}

double BinLog::dudg(const double g) const {
	return 1. / (g * log(10.));
}
