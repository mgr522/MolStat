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
#include <cmath>
#include <valarray>
#include <iostream>
#include <fstream>

#include <general/string_tools.h>
#include <general/random_distributions/rng.h>
#include <general/histogram_tools/counterindex.h>
#include <general/histogram_tools/histogram.h>

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

	// open the output file
	ofstream histout(parser.outputFileName());
	if(!histout)
	{
		cout << "FATAL ERROR: Unable to open \"" << parser.outputFileName() <<
			"\" for output." << endl;
		return 0;
	}

	// print the simulator information
	parser.printState(cout);

	// initialize the GSL random number generator
	gsl_rng_env_setup();
	molstat::gsl_rng_ptr r(gsl_rng_alloc(gsl_rng_default), &gsl_rng_free);
	//gsl_rng_set(r.get(), 0xFEEDFACE); // use this line for debugging
	gsl_rng_set(r.get(), time(nullptr));

	// create the histogram object
	// first need the bin styles to determine the dimensionality
	vector<shared_ptr<const molstat::BinStyle>> bstyles(0);
	{
		const auto nonconst = parser.getBinStyles();
		bstyles.resize(nonconst.size());
		for(size_t j = 0; j < nonconst.size(); ++j)
			bstyles[j] = nonconst[j];
	} // this was necessary to add const to the pointer

	molstat::Histogram hist(bstyles.size());

	// Get the requested number of samples
	for (size_t j = 0; j < ntrials; ++j)
	{
		// add the data to the list
		hist.add_data(sim->simulate(r));
	}

	// make the histogram
	hist.bin_data(bstyles);
	// DO ERROR CHECKING ----------------------------------------------------

	// go through the bins
	molstat::CounterIndex ci { hist.begin() };

	for(molstat::CounterIndex ci{ hist.begin() }; !ci.at_end(); ++ci)
	{
		const valarray<double> coords = hist.getCoordinates(ci);
		for(size_t j = 0; j < coords.size(); ++j)
			histout << coords[j] << ' ';
		histout << hist.getBinCount(ci) << endl;
	}

	histout.close();

	return 0;
}