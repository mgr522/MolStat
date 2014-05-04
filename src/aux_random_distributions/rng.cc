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
#include "../string_tools.h"
#include "constant.h"
#include "uniform.h"
#include "normal.h"

using namespace std;

shared_ptr<RandomDistribution> distribution_from_tokens(
	const std::vector<std::string> &tokens) {

	shared_ptr<RandomDistribution> ret;

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

		try {
			val = stod(tokens[1]);
		}
		catch(const invalid_argument &e) {
			throw invalid_argument("Constant Distribution: Unable to convert " \
				"the first argument to a double.");
		}

		ret = make_shared<ConstantDistribution>(val);
	}
	else if(type == "uniform") {
		double lower, upper;

		// tokens[1] and tokens[2] are the lower and upper bounds of the range,
		// respectively
		if(tokens.size() < 3)
			throw invalid_argument("Invalid uniform distribution. Use\n" \
				"   uniform lower upper\n" \
				"where lower and upper are the bounds, respectively.");

		try {
			lower = stod(tokens[1]);
			upper = stod(tokens[2]);
		}
		catch(const invalid_argument &e) {
			throw invalid_argument("Uniform Distribution: Unable to convert " \
				"the first and/or second argument to a double.");
		}

		if(lower >= upper)
			throw invalid_argument("Uniform Distribution: The lower bound must " \
				"be lower than the upper bound.");

		ret = make_shared<UniformDistribution>(lower, upper);
	}
	else if(type == "normal" || type == "gaussian") {
		double mean, stdev;

		// tokens[1] and tokens[2] are the mean and standard deviation of the
		// distribution, respectively
		if(tokens.size() < 3)
			throw invalid_argument("Invalid normal distribution. Use\n" \
				"   normal mean standard-deviation\n");

		try {
			mean = stod(tokens[1]);
			stdev = stod(tokens[2]);
		}
		catch(const invalid_argument &e) {
			throw invalid_argument("Normal Distribution: Unable to convert " \
				"the first and/or second argument to a double.");
		}

		if(stdev <= 0.)
			throw invalid_argument("Normal Distribution: The standard deviation" \
				"must be positive.");

		ret = make_shared<NormalDistribution>(mean, stdev);
	}
	else
		throw invalid_argument("Unrecognized probability distribution.\n" \
			"Possible options are:\n" \
			"   Constant - Specify a value.\n" \
			"   Uniform - Uniform distribution.\n" \
			"   Normal - Normal (Gaussian) distribution.\n" \
			"   Gaussian - Normal (Gaussian) distribution.\n");

	return ret;
}
