/* This file is a part of MolStat, which is distributed under the Creative
   Commons Attribution-NonCommercial 4.0 International Public License.
   MolStat (c) 2014, Northwestern University. */
/**
 * \file rng.cc
 * \brief Implementation of the function for generating a random number
 *    distribution from a tokenized line of input.
 *
 * \author Matthew G.\ Reuter
 * \date May 2014
 */

#include "rng.h"
#include <gsl/gsl_randist.h>
#include <general/string_tools.h>
#include "constant.h"
#include "uniform.h"
#include "normal.h"
#include "lognormal.h"
#include "gamma.h"

namespace molstat {

using namespace std;

std::unique_ptr<RandomDistribution> RandomDistributionFactory(
	const std::vector<std::string> &tokens) {

	unique_ptr<RandomDistribution> ret;

	if(tokens.size() < 1)
		throw invalid_argument("Empty line.");

	// the first token is the type of distribution to form
	string type = tokens[0];
	make_lower(type);

	if(type == "constant") {
		double val;

		// tokens[1] is the value to return; cast it to a double
		if(tokens.size() < 2)
			throw invalid_argument("Invalid constant distribution. Use\n" \
				"   constant value\n" \
				"where value is the value to be returned.");

		val = string_to_double(tokens[1]);

		ret = unique_ptr<RandomDistribution>(new ConstantDistribution(val));
	}
	else if(type == "uniform") {
		double lower, upper;

		// tokens[1] and tokens[2] are the lower and upper bounds of the range,
		// respectively
		if(tokens.size() < 3)
			throw invalid_argument("Invalid uniform distribution. Use\n" \
				"   uniform lower upper\n" \
				"where lower and upper are the bounds, respectively.");

		lower = string_to_double(tokens[1]);
		upper = string_to_double(tokens[2]);

		if(lower >= upper)
			throw invalid_argument("Uniform Distribution: The lower bound must " \
				"be lower than the upper bound.");

		ret = unique_ptr<RandomDistribution>(
			new UniformDistribution(lower, upper));
	}
	else if(type == "normal" || type == "gaussian") {
		double mean, stdev;

		// tokens[1] and tokens[2] are the mean and standard deviation of the
		// distribution, respectively
		if(tokens.size() < 3)
			throw invalid_argument("Invalid normal distribution. Use\n" \
				"   normal mean standard-deviation");

		mean = string_to_double(tokens[1]);
		stdev = string_to_double(tokens[2]);

		if(stdev <= 0.)
			throw invalid_argument("Normal Distribution: The standard deviation" \
				" must be positive.");

		ret = unique_ptr<RandomDistribution>(
			new NormalDistribution(mean, stdev));
	}
	else if(type == "lognormal") {
		double zeta, sigma;

		// tokens[1] and tokens[2] are the mean and standard deviation (in log-
		// space), respectively
		if(tokens.size() < 3)
			throw invalid_argument("Invalid lognormal distribution. Use\n" \
				"   lognormal zeta sigma");

		zeta = string_to_double(tokens[1]);
		sigma = string_to_double(tokens[2]);

		if(sigma <= 0.)
			throw invalid_argument("Lognormal Distribution: The standard " \
				"deviation (sigma) must be positive.");

		ret = unique_ptr<RandomDistribution>(
			new LognormalDistribution(zeta, sigma));
	}
	else if(type == "gamma") {
		double shape, scale;

		// tokens[1] and tokens[2] are the shape and scale factors, respectively
		if(tokens.size() < 3)
			throw invalid_argument("Invalid gamma distribution. Use\n" \
				"   gamma shape scale");

		shape = string_to_double(tokens[1]);
		scale = string_to_double(tokens[2]);

		if(shape <= 0. || scale <= 0.)
			throw invalid_argument("Gamma Distribution: The shape and scale " \
				"factors must be positive.");

		ret = unique_ptr<RandomDistribution>(
			new GammaDistribution(shape, scale));
	}
	else
		throw invalid_argument("Unrecognized probability distribution.\n" \
			"Possible options are:\n" \
			"   Constant - Specify a value.\n" \
			"   Uniform - Uniform distribution.\n" \
			"   Normal - Normal (Gaussian) distribution.\n" \
			"   Gaussian - Normal (Gaussian) distribution.\n" \
			"   Lognormal - Lognormal distribution.\n" \
			"   Gamma - Gamma distribution.\n");

	return ret;
}

} // namespace molstat