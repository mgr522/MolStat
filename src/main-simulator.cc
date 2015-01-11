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
#include <general/histogram_tools/bin_linear.h>
#include <general/simulator_tools/simulator_exceptions.h>

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

	// initialize the random number engine
	//molstat::Engine engine{ 0xFEEDFACE }; // use this line for debugging
	molstat::Engine engine{ static_cast<unsigned int>(time(nullptr)) };

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
	// count the number of trials that don't emit the observable
	size_t no_obs { 0 };
	for (size_t j = 0; j < ntrials; ++j)
	{
		try
		{
			// add the data to the list
			hist.add_data(sim->simulate(engine));
		}
		catch(const molstat::NoObservableProduced &e)
		{
			// one of the observables was not emitted for the randomly generated
			// parameters
			++no_obs;
		}
	}

	// print out the number of trials that did not produce an observable
	cout << '\n' << no_obs << " of the " << ntrials << " trials (" <<
		(100. * no_obs / ntrials) << "%) did not produce an observable." << endl;

	// make the histogram
	// if we encounter a bad dimension -- specifically, one where there is no
	// range of data (all trials yield the same value) and more than one bin
	// is specified -- override the binstyle for that dimension and try again
	bool binned { false };
	do
	{
		try
		{
			hist.bin_data(bstyles);
			binned = true;
		}
		catch(const size_t &bad_dim)
		{
			// one dimension specified multiple bins and does not have a range of
			// values
			cout << "Empty data range in dimension " << bad_dim << "; however, " \
				"more than 1 bin was requested.\nOnly using 1 bin." << endl;
			bstyles[bad_dim] = make_shared<const molstat::BinLinear>(1);
		}
	} while(!binned);

	// go through the bins
	molstat::CounterIndex ci { hist.begin() };

	for(molstat::CounterIndex ci{ hist.begin() }; !ci.at_end(); ++ci)
	{
		const valarray<double> coords = hist.getCoordinates(ci);
		for(size_t j = 0; j < coords.size(); ++j)
			histout << coords[j] << ' ';
		histout << hist.getBinCount(ci) << endl;
	}

	// close the output stream
	histout.close();

	return 0;
}
