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
#include "electron_transport/simulator_models/transport_simulate_models.h"

using namespace std;

/**
 * \internal
 * \brief Enum for the type of histogram (1D or 2D).
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
 * \todo Add a mechanism to specify the binning style of each axis.
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

	// list of parameters with specified distributions
	map<string, shared_ptr<RandomDistribution>> parameters;

	// the model
	shared_ptr<SimulateModel> model;
	// the list of models:
	// stored as a map of string (model name) to a function that instantiates
	// such a model, given a map of random number distributions.
	map<string, SimulateModelInstantiator> models;

	// function for generating our observable
	Observable<2> observable;
	// the list of available observables:
	// stored as a map of string (observable name) to a function that takes
	// a SimulateModel and produces the observable function.
	// This extra function layer is needed to dynamically typecheck the model,
	// making sure the model/observable pair is valid.
	map<string, function< Observable<2>(shared_ptr<SimulateModel>)> >
		observables;

	// I/O variables
	string line, modelname, obsname, name;
	vector<string> tokens;
	shared_ptr<BinStyle> bstyle;

	// variables for setting up the histograms / storing the random data
	array<double, 2>
		mins{{numeric_limits<double>::max(), numeric_limits<double>::max()}},
		maxs{{-numeric_limits<double>::max(), -numeric_limits<double>::max()}};
	forward_list<array<double, 2>> data;

	// load models and observables
	// SimulateModelAdd calls appear here
	load_transport_models(models);
	// ObservableCheck calls appear here
	load_transport_observables(observables);

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
	modelname = tokens[0];

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
	// store the observable name for later.
	obsname = tokens[0];

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
	// first token: number of bins to use
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

	// try to instantiate the model
	try {
		// models.at() returns a function for instantiating the model
		model = models.at(to_lower(modelname))(parameters);
	}
	catch(const out_of_range &e) {
		// model not found
		fprintf(stderr, "Error: model \"%s\" not found.\n", modelname.c_str());
		return 0;
	}
	catch(const runtime_error &e) {
		// at least one required parameter was unspecified
		fprintf(stderr, "Error: a distribution for \"%s\" must be specified.\n",
			e.what());
		return 0;
	}

	// now verify the observable
	try {
		// observables.at() returns a function to check the model/observable
		// combination; the second function actually does the check.
		observable = observables.at(to_lower(obsname))(model);
	}
	catch(const out_of_range &e) {
		// observable not found
		fprintf(stderr, "Error: observable \"%s\" not found.\n",
			obsname.c_str());
		return 0;
	}
	catch(const runtime_error &e) {
		fprintf(stderr, "Error: model \"%s\" and observable \"%s\" are " \
			"incompatible.\n", modelname.c_str(), obsname.c_str());
		return 0;
	}

	// Get the requested number of samples
	for (i = 0; i < ntrials; ++i) {
		array<double, 2> datum;

		datum = observable(r);

		// check the limits
		if(datum[0] < mins[0])
			mins[0] = datum[0];
		if(datum[0] > maxs[0])
			maxs[0] = datum[0];
		if(datum[1] < mins[1])
			mins[1] = datum[1];
		if(datum[1] > maxs[1])
			maxs[1] = datum[1];

		// add the data to the list
		data.emplace_front(datum);
	}

	// histogram variables
	HistogramType htype;
	shared_ptr<Histogram1D> hist1;
	shared_ptr<Histogram2D> hist2;
	function<double(const array<double, 2> &)> mask1;
	array<double, 2> bounds1;

	// setup the histogram type; what type are we making?
	htype = HistogramType::TwoD; // base guess

	// check the range of the first variable
	if(mins[0] == maxs[0] && mins[1] == maxs[1]) {
		fprintf(stderr, "Error: no data to bin into histograms. All simulated " \
			"data have the value\n       {%f, %f}\n", mins[0], mins[1]);
		return 0;
	}
	else if(mins[0] == maxs[0]) {
		// process only the second variable
		htype = HistogramType::OneD;
		mask1 = [] (const array<double, 2> &a) -> double {
			return a[1];
		};
		bounds1[0] = mins[1];
		bounds1[1] = maxs[1];
	}
	else if(mins[1] == maxs[1]) {
		// process only the second variable
		htype = HistogramType::OneD;
		mask1 = [] (const array<double, 2> &a) -> double {
			return a[0];
		};
		bounds1[0] = mins[0];
		bounds1[1] = maxs[0];
	}
	else {
		// add a pinch to the max since the GSL histogrammer is exclusive on the
		// upper bounds of binning ranges
		maxs[0] += 1.e-6;
		maxs[1] += 1.e-6;
	}

	// set up the histogram, populate it, and print it
	if(htype == HistogramType::OneD) {
		// note the similar pinch added to bounds[1] for the above reason
		hist1 = make_shared<Histogram1D>(nbin, bounds1[0], bounds1[1] + 1.e-6,
			bstyle);

		while(!data.empty()) {
			hist1->add_data(mask1(data.front()));
			data.pop_front();
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

		while(!data.empty()) {
			hist2->add_data(data.front());
			data.pop_front();
		}

		for(Histogram2D::const_iterator iter = hist2->begin();
			iter != hist2->end();
			++iter) {

			printf("%.6f %.6f %.6f\n", iter.get_variable()[0],
				iter.get_variable()[1], iter.get_bin_count());
		}
	}

	return 0;
}
