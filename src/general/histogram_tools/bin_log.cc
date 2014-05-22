/* This file is a part of MolStat, which is distributed under the Creative
   Commons Attribution-NonCommercial 4.0 International Public License.
   MolStat (c) 2014, Northwestern University. */
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
