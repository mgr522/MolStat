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
 * \date October 2014
 */

#include <stdexcept>
#include <functional>
#include <cmath>
#include <forward_list>
#include <valarray>
#include <limits>

#include <general/string_tools.h>
#include <general/random_distributions/rng.h>
#include <general/histogram_tools/histogram1d.h>
#include <general/histogram_tools/histogram2d.h>

#include "main-simulator.h"

using namespace std;

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
int main(int argc, char **argv)
{
	// name of the histogram output file with a default name
	string histfilename{"histogram.dat"};

	// process the input deck
	SimulatorInputParse parser;
	try
	{
		parser.readInput(cin, cout);
	}
	catch(const runtime_error &e)
	{
		cout << "FATAL ERROR: " << e.what() << endl;
		return 0;
	}

	// for debugging purposes, we may want to print the state here
	// parser.printState(cout);

	// readInput did some checking (syntax of input file, etc.), but is
	// incomplete... no verifying model names are correct, all observables
	// and binning styles are specified, etc.

	// need to finish processing input
	// check the number of trials
	const size_t ntrials{ parser.numTrials() };
	if(ntrials == 0)
	{
		cout << "FATAL ERROR: There must be at least one trial." << endl;
		return 0;
	}

	// create the simulator
	// this will make sure model names are good, all distributions are
	// specified, etc.
	unique_ptr<molstat::Simulator> sim{ nullptr };
	try
	{
		sim = parser.createSimulator(cout);
	}
	catch(const exception &e)
	{
		cout << "FATAL ERROR: " << e.what() << endl;
		return 0;
	}

	// print the simulator information
	parser.printState(cout);

#if 0
	// initialize the GSL random number generator
	gsl_rng_env_setup();
	molstat::gsl_rng_ptr r(gsl_rng_alloc(gsl_rng_default), &gsl_rng_free);
	//gsl_rng_set(r.get(), 0xFEEDFACE); // use this line for debugging
	gsl_rng_set(r.get(), time(nullptr));

	// simulation variables
	size_t nbin, i;
	shared_ptr<molstat::BinStyle> bstyle;

	// variables for setting up the histograms / storing the random data
	array<double, 2>
		mins{{numeric_limits<double>::max(), numeric_limits<double>::max()}},
		maxs{{-numeric_limits<double>::max(), -numeric_limits<double>::max()}};
	forward_list<array<double, 2>> data;

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
	unique_ptr<molstat::Histogram1D> hist1;
	unique_ptr<molstat::Histogram2D> hist2;
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
		hist1.reset(new molstat::Histogram1D(nbin, bounds1[0],
			bounds1[1] + 1.e-6, bstyle));

		while(!data.empty()) {
			hist1->add_data(mask1(data.front()));
			data.pop_front();
		}

		for(molstat::Histogram1D::const_iterator iter = hist1->begin();
			iter != hist1->end();
			++iter) {

			printf("%.6f %.6f\n", iter.get_variable()[0], iter.get_bin_count());
		}
	}
	else if(htype == HistogramType::TwoD) {
		hist2.reset(new molstat::Histogram2D(array<size_t, 2>{{nbin, nbin}},
			mins, maxs, bstyle));

		while(!data.empty()) {
			hist2->add_data(data.front());
			data.pop_front();
		}

		for(molstat::Histogram2D::const_iterator iter = hist2->begin();
			iter != hist2->end();
			++iter) {

			printf("%.6f %.6f %.6f\n", iter.get_variable()[0],
				iter.get_variable()[1], iter.get_bin_count());
		}
	}

#endif
	return 0;
}