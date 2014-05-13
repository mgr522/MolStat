/**
 * \file bin_log.cc
 * \brief Implements logarithmic binning.
 *
 * \author Matthew G.\ Reuter
 * \date May 2014
 */

#include "bin_log.h"
#include <cmath>

BinLog::BinLog(const double b_)
	: b(b_) {
}

double BinLog::gmask(const double g) const {
	return log(g) / log(b);
}

double BinLog::invgmask(const double u) const {
	return pow(b, u);
}

double BinLog::dudg(const double g) const {
	return 1. / (g * log(b));
}
