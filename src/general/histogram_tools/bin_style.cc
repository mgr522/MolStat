/* This file is a part of MolStat, which is distributed under the Creative
   Commons Attribution-NonCommercial 4.0 International Public License.
   MolStat (c) 2014, Northwestern University. */
/**
 * \file bin_style.cc
 * \brief Implements a function to create a BinStyle from a vector of tokens.
 *
 * \author Matthew G.\ Reuter
 * \date June 2014
 */

#include "bin_style.h"
#include "bin_linear.h"
#include "bin_log.h"
#include "../string_tools.h"

using namespace std;

namespace molstat {

std::unique_ptr<BinStyle> BinStyleFactory(
	const std::vector<std::string> &tokens) {

	unique_ptr<BinStyle> ret;

	if(tokens.size() < 1)
		throw invalid_argument("Empty line.");

	// the first token is the name of the binning style
	string name = tokens[0];
	make_lower(name);

	if(name == "linear") {
		ret.reset(new BinLinear());
	}
	else if(name == "log") {
		// need to read the base, if available. If not, use 10.
		double b;

		if(tokens.size() >= 2) {
			b = string_to_double(tokens[1]);
			if(b <= 0.)
				throw invalid_argument("The logarithm base must be positive.");
		}
		else
			b = 10.;

		ret.reset(new BinLog(b));
	}
	else
		throw invalid_argument("Unrecognized binning style.\n" \
			"Possible options are:\n" \
			"   Linear - Linear binning.\n" \
			"   Log - Logarithmic binning (base defaults to 10).\n");

	return ret;
}

} // namespace molstat
