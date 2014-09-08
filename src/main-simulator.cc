/* This file is a part of MolStat, which is distributed under the Creative
   Commons Attribution-NonCommercial 4.0 International Public License.
   MolStat (c) 2014, Northwestern University. */
/**
 * \file main-simulator.cc
 * \brief Main function for simulating histograms of single-molecule
 *    properties.
 *
 * This code reads in various input parameters from standard in, and uses these
 * parameters to simulate data. The data is then binned into 1D or 2D
 * histograms, depending on the desired observable.
 *
 * \author Matthew G.\ Reuter
 * \date September 2014
 */

#include <memory>
#include <cstdio>
#include <stdexcept>
#include <string>
#include <vector>
#include <functional>
#include <map>
#include <cmath>
#include <forward_list>
#include <array>
#include <limits>

#include <general/string_tools.h>
#include <general/random_distributions/rng.h>
#include <general/histogram_tools/histogram1d.h>
#include <general/histogram_tools/histogram2d.h>
#include <general/histogram_tools/bin_linear.h>
#include <general/histogram_tools/bin_log.h>
//#include "simulator_models/conductance_model.h"

using namespace std;

/**
 * \internal
 * \brief Enum for the type of histogram (1D or 2D).
 *
 * A `ZeroBias` CalculationType is always 1D. `Static` and `Differential` are
 * usually 2D, but can be 1D if the voltage's random distribution is a
 * ConstantDistribution.
 * \endinternal
 */
enum class HistogramType {
	OneD, // 1D Histogram
	TwoD // 2D Histogram
};

/**
 * \internal
 * \brief Main function for simulating a histogram.
 *
 * Parses the input parameters and outputs randomly generated data
 * for the desired observable.
 *
 * \param[in] argc The number of command-line arguments.
 * \param[in] argv The command-line arguments.
 * \return Exit status; 0 for normal.
 * \endinternal
 */
int main(int argc, char **argv) {
	// initialize the GSL random number generator
	gsl_rng_env_setup();
	shared_ptr<gsl_rng> r(gsl_rng_alloc(gsl_rng_default), &gsl_rng_free);
	//gsl_rng_set(r.get(), 0xFEEDFACE); // use this line for debugging
	gsl_rng_set(r.get(), time(nullptr));

	// simulation variables
	size_t ntrials, nbin, i;

	// model variables
	map<string, shared_ptr<RandomDistribution>> parameters;
	//shared_ptr<ConductanceModel> model;
	shared_ptr<RandomDistribution> dist_V, dist_eta;

	// I/O variables
	string line, modeltype, name;
	vector<string> tokens;

	// histogram variables
	shared_ptr<BinStyle> bstyle;
	//CalculationType ctype;
	HistogramType htype;
	shared_ptr<Histogram1D> hist1;
	shared_ptr<Histogram2D> hist2;

	// variables for setting up the histograms / storing the random data
	array<double, 2>
		mins{{numeric_limits<double>::max(), numeric_limits<double>::max()}},
		maxs{{-numeric_limits<double>::max(), -numeric_limits<double>::max()}};
	forward_list<array<double, 2>> gdata;

#if 0
	// setup the simulation -- read in parameters from stdin
	// Line 1: One token specifying the model to use
	try {
		line = getline(stdin);
	}
	catch(const runtime_error &e) {
		fprintf(stderr, "Error: %s.\n", e.what());
		return 0;
	}

	tokenize(line, tokens);
	if(tokens.size() < 1) {
		fprintf(stderr, "Error: model name expected in line 1.\n");
		return 0;
	}

	// store the model name for later (need to continue reading data before
	// invoking the constructor.
	modeltype = tokens[0];
	make_lower(modeltype);

	// Line 2: Observable
	try {
		line = getline(stdin);
	}
	catch(const runtime_error &e) {
		fprintf(stderr, "Error: %s.\n", e.what());
		return 0;
	}
	tokenize(line, tokens);
	if(tokens.size() < 1) {
		fprintf(stderr, "Error: conductance type expected in line 2.\n");
		return 0;
	}
	make_lower(tokens[0]);
	if(tokens[0] == "static") {
		ctype = CalculationType::Static;
		htype = HistogramType::TwoD;
	}
	else if(tokens[0] == "differential") {
		ctype = CalculationType::Differential;
		htype = HistogramType::TwoD;
	}
	else if(tokens[0] == "zerobias") {
		ctype = CalculationType::ZeroBias;
		htype = HistogramType::OneD;
	}
	else {
		fprintf(stderr, "Error: Unrecognized conductance type in line 2.\n   " \
			"It must be \"Static\", \"Differential\", or \"ZeroBias\".\n");
		return 0;
	}

	// Line 3: number of trials
	try {
		line = getline(stdin);
	}
	catch(const runtime_error &e) {
		fprintf(stderr, "Error: %s.\n", e.what());
		return 0;
	}
	tokenize(line, tokens);
	if(tokens.size() < 1) {
		fprintf(stderr, "Error: number of trials expected in line 3.\n");
		return 0;
	}
	try {
		ntrials = stoul(tokens[0]);
	}
	catch(const invalid_argument &e) {
		fprintf(stderr, "Error: unrecognizable number '%s'.\n",
			tokens[0].c_str());
		return 0;
	}

	// Line 4: Binning information
	try {
		line = getline(stdin);
	}
	catch(const runtime_error &e) {
		fprintf(stderr, "Error: %s.\n", e.what());
		return 0;
	}
	tokenize(line, tokens);
	if(tokens.size() < 2) {
		fprintf(stderr, "Error: binning information expected in line 4.\n" \
			"Number of bins followed by the type of bins.\n");
		return 0;
	}

	// first: number of bins to use
	try {
		nbin = stoul(tokens[0]);
	}
	catch(const invalid_argument &e) {
		fprintf(stderr, "Error: unrecognizable number '%s'.\n",
			tokens[0].c_str());
		return 0;
	}
	// remove the number of bins so we can get the binning style
	tokens.erase(tokens.begin());
	make_lower(tokens[0]);
	try {
		bstyle = get_bin_style(tokens);
	}
	catch(const invalid_argument &e) {
		fprintf(stderr, "Error: unknown binning style.\n");
		return 0;
	}
#endif

	// all subsequent lines specify random number distributions
	// EOF is flagged by a runtime_error in the getline function
	try {
		while(true) {
			line = getline(stdin);

			tokenize(line, tokens);
			if(tokens.size() > 0) {
				name = tokens[0];
				make_lower(name);

				// need to remove the front entry for the distribution generator
				tokens.erase(tokens.begin());

				// create a random distribution from the contents of the line
				try {
					// create the separate line; if distribution_from_tokens throws
					// an exception parameters would still create an entry for name
					shared_ptr<RandomDistribution> rd =
						distribution_from_tokens(tokens);
					parameters[name] = rd;
				}
				catch(const invalid_argument &e) {
					fprintf(stderr, "Error: unable to form a random number " \
						"distribution from:\n   %s\n%s\n", line.c_str(), e.what());
				}
			}
		}
	}
	catch(const runtime_error &e) {
		// this just means we hit EOF -- stop trying to read more
	}

#if 0
	// set the distributions required by the model
	try {
		model = make_model(modeltype, parameters);
	}
	catch(const invalid_argument &e) {
		fprintf(stderr, "Error initializing the model: %s\n", e.what());
		return 0;
	}

	// set up other variables for the simulation, including setting other
	// distributions
	if(ctype == CalculationType::Static ||
		ctype == CalculationType::Differential) {

		// make sure there are distributions for V and eta
		try {
			dist_eta = parameters.at("eta");
		}
		catch(const out_of_range &e) {
			fprintf(stderr, "Error: a distribution for \"eta\" must be " \
				"specified.\n");
			return 0;
		}
		try {
			dist_V = parameters.at("v");
		}
		catch(const out_of_range &e) {
			fprintf(stderr, "Error: a distribution for \"V\" must be specified." \
				"\n");
			return 0;
		}

		// check if V is a ConstantDistribution. If it is, this is really a 1D
		// histogram
		if(dynamic_pointer_cast<ConstantDistribution>(dist_V) != nullptr)
			htype = HistogramType::OneD;
	}
	else if(ctype == CalculationType::ZeroBias) {
		// nothing required
	}
	else {
		// should never be here
		fprintf(stderr, "An unknown error has occured.\n");
		return 0;
	}

	// Get the requested number of voltage/conductance data points
	for (i = 0; i < ntrials; ++i) {
		double V, GV;

		if(ctype == CalculationType::Static ||
			ctype == CalculationType::Differential) {

			V = dist_V->sample(r);

			if(ctype == CalculationType::Static)
				GV = model->static_conductance(r, EF, dist_eta->sample(r), V);
			else
				GV = model->diff_conductance(r, EF, dist_eta->sample(r), V);
		}
		else if(ctype == CalculationType::ZeroBias) {
			V = 0.;
			GV = model->zero_bias_conductance(r, EF);
		}

		// check the limits
		if(V < mins[0])
			mins[0] = V;
		if(V > maxs[0])
			maxs[0] = V;
		if(GV < mins[1])
			mins[1] = GV;
		if(GV > maxs[1])
			maxs[1] = GV;

		// add the data to the list
		gdata.emplace_front(array<double, 2>{{V, GV}});
	}

	// set up the histogram, populate it, and print it
	if(htype == HistogramType::OneD) {
		hist1 = make_shared<Histogram1D>(nbin, mins[1], maxs[1], bstyle);

		while(!gdata.empty()) {
			hist1->add_data(gdata.front()[1]);
			gdata.pop_front();
		}

		for(Histogram1D::const_iterator iter = hist1->begin();
			iter != hist1->end();
			++iter) {

			printf("%.6f %.6f\n", iter.get_variable()[0], iter.get_bin_count());
		}
	}
	else if(htype == HistogramType::TwoD) {
		hist2 = make_shared<Histogram2D>(array<size_t, 2>{{nbin, nbin}},
			mins, maxs, bstyle);

		while(!gdata.empty()) {
			hist2->add_data(gdata.front());
			gdata.pop_front();
		}

		for(Histogram2D::const_iterator iter = hist2->begin();
			iter != hist2->end();
			++iter) {

			printf("%.6f %.6f %.6f\n", iter.get_variable()[0],
				iter.get_variable()[1], iter.get_bin_count());
		}
	}
#endif

	return 0;
}
