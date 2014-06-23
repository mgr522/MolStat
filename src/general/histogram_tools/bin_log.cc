/* This file is a part of MolStat, which is distributed under the Creative
   Commons Attribution-NonCommercial 4.0 International Public License.
   MolStat (c) 2014, Northwestern University. */
/**
 * \internal
 * \file bin_log.cc
 * \brief Implements logarithmic binning.
 *
 * \author Matthew G.\ Reuter
 * \date June 2014
 * \endinternal
 */

#include "bin_log.h"
#include <cmath>

BinLog::BinLog(const double b_)
	: b(b_) {
}

double BinLog::mask(const double x) const {
	return log(x) / log(b);
}

double BinLog::invmask(const double u) const {
	return pow(b, u);
}

double BinLog::dmaskdx(const double x) const {
	return 1. / (x * log(b));
}
