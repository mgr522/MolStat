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
#include <general/simulator_tools/simulator_exceptions.h>
#include <electron_transport/simulator_models/transport_simulate_module.h>

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

	// process the input deck and get our simulator
	unique_ptr<molstat::Simulator> sim;
	try {
		sim = processInput(cin, histfilename);
	}
	catch(const exception &e) {
		cout << "FATAL ERROR: " << e.what() << endl;
		return 0;
	}

#if 0
	// initialize the GSL random number generator
	gsl_rng_env_setup();
	molstat::gsl_rng_ptr r(gsl_rng_alloc(gsl_rng_default), &gsl_rng_free);
	//gsl_rng_set(r.get(), 0xFEEDFACE); // use this line for debugging
	gsl_rng_set(r.get(), time(nullptr));

	// simulation variables
	size_t ntrials, nbin, i;

	// list of parameters with specified distributions
	map<string, shared_ptr<molstat::RandomDistribution>> parameters;

	// the model
	shared_ptr<molstat::Simulator<2>> model;

	// function for generating our observable
	molstat::Observable<2> observable;

	// I/O variables
	string line, modelname, obsname, name;
	vector<string> tokens;
	shared_ptr<molstat::BinStyle> bstyle;

	// variables for setting up the histograms / storing the random data
	array<double, 2>
		mins{{numeric_limits<double>::max(), numeric_limits<double>::max()}},
		maxs{{-numeric_limits<double>::max(), -numeric_limits<double>::max()}};
	forward_list<array<double, 2>> data;

	// setup the simulation -- read in parameters from stdin
	// Line 1: One token specifying the model to use
	try {
		line = molstat::getline(stdin);
	}
	catch(const runtime_error &e) {
		fprintf(stderr, "Error: %s.\n", e.what());
		return 0;
	}

	molstat::tokenize(line, tokens);
	if(tokens.size() < 1) {
		fprintf(stderr, "Error: model name expected in line 1.\n");
		return 0;
	}
	// store the model name for later (need to continue reading data before
	// invoking the constructor.
	modelname = tokens[0];

	// Line 2: Observable
	try {
		line = molstat::getline(stdin);
	}
	catch(const runtime_error &e) {
		fprintf(stderr, "Error: %s.\n", e.what());
		return 0;
	}
	molstat::tokenize(line, tokens);
	if(tokens.size() < 1) {
		fprintf(stderr, "Error: conductance type expected in line 2.\n");
		return 0;
	}
	// store the observable name for later.
	obsname = tokens[0];

	// Line 3: number of trials
	try {
		line = molstat::getline(stdin);
	}
	catch(const runtime_error &e) {
		fprintf(stderr, "Error: %s.\n", e.what());
		return 0;
	}
	molstat::tokenize(line, tokens);
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
		line = molstat::getline(stdin);
	}
	catch(const runtime_error &e) {
		fprintf(stderr, "Error: %s.\n", e.what());
		return 0;
	}
	molstat::tokenize(line, tokens);
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
	molstat::make_lower(tokens[0]);
	try {
		bstyle = molstat::BinStyleFactory(tokens);
	}
	catch(const invalid_argument &e) {
		fprintf(stderr, "Error: unknown binning style.\n");
		return 0;
	}

	// all subsequent lines specify random number distributions
	// EOF is flagged by a runtime_error in the getline function
	try {
		while(true) {
			line = molstat::getline(stdin);

			molstat::tokenize(line, tokens);
			if(tokens.size() > 0) {
				name = tokens[0];
				molstat::make_lower(name);

				// need to remove the front entry for the distribution generator
				tokens.erase(tokens.begin());

				// create a random distribution from the contents of the line
				try {
					// create the separate line; if distribution_from_tokens throws
					// an exception parameters would still create an entry for name
					shared_ptr<molstat::RandomDistribution> rd =
						molstat::RandomDistributionFactory(tokens);
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
		model = models.at(molstat::to_lower(modelname))(parameters);
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
		observable = observables.at(molstat::to_lower(obsname))(model);
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

// ************************************************************************
// BEGIN INPUT AND COMMAND PROCESSING FUNCTIONS
// ************************************************************************
std::unique_ptr<molstat::Simulator> processInput(std::istream &input,
	std::string &ofname)
{
	// the list of models:
	// stored as a map of string (model name) to a function that produces a
	// factory for that model.
	map<string, molstat::SimulateModelFactoryFunction> models;

	// the list of available observables:
	// stored as a map of string (observable name) to the index of the
	// observable.
	map<string, molstat::ObservableIndex> observables;

	// load models and observables
	// the general syntax for loading models is
	// models.emplace( to_lower(name),
	//                 molstat::GetSimulateModelFactory<model_type>() );
	//
	// and likewise for observables,
	// observables.emplace( to_lower(name),
	//                      GetObservableIndex<observable_type>() );
	molstat::transport::load_models(models);
	molstat::transport::load_observables(observables);

	// a pointer to the model that will get fed to the simulator
	shared_ptr<molstat::SimulateModel> model{ nullptr };

	// a map of observables, from number to index
	map<size_t, molstat::ObservableIndex> used_observables;

	// process the input deck
	// input evaluates to "false" when we hit EOF
	while(input)
	{
		string line;
		getline(input, line);

		// tokenize the string
		queue<string> tokens = molstat::tokenize(line);
		if(tokens.size() == 0) // empty line
			continue;

		// the first token is the command name, pop it off and then process the
		// rest of the tokens
		string command { molstat::to_lower(tokens.front()) };
		tokens.pop();

		// go through the supported commands
		if(command == "model")
		{
			try
			{
				model = processModel(move(tokens), input, models);
			}
			catch(const runtime_error &e) {
				cout << "Error processing \"" + line + "\":\n   " << e.what() <<
					endl;
			}
		}
		else if(command == "observable" || command == "observable_x")
		{
			try
			{
				used_observables.emplace(
					0, processObservable(move(tokens), observables));
			}
			catch(const runtime_error &e)
			{
				cout << "Error processing \"" + line + "\":\n   " << e.what() <<
					endl;
			}
		}
		else if(command == "observable_y")
		{
			try
			{
				used_observables.emplace(
					1, processObservable(move(tokens), observables));
			}
			catch(const runtime_error &e)
			{
				cout << "Error processing \"" + line + "\":\n   " << e.what() <<
					endl;
			}
		}
		else
		{
			cout << "Unrecognized command: \"" << command << "\"." << endl;
		}
	}

	// construct the Simulator
	unique_ptr<molstat::Simulator> sim{ new molstat::Simulator(model) };

	// load the observables
	for(const auto obs : used_observables)
	{
		try
		{
			sim->setObservable(obs.first, obs.second);
		}
		catch(const out_of_range &e)
		{
			throw runtime_error(
				"Non-contiguous observable numbering is not allowed.");
		}
		catch(const molstat::IncompatibleObservable &e)
		{
			throw runtime_error("The model is incompatible with observable "
				+ to_string(obs.first) + '.');
		}
	}

	return sim;
}

std::shared_ptr<molstat::SimulateModel> processModel(
	std::queue<std::string> &&tokens,
	std::istream &input,
	const std::map<std::string,
	               molstat::SimulateModelFactoryFunction> &models)
{
	if(tokens.size() == 0)
		throw runtime_error("Model name not specified.");

	// get the factory
	molstat::SimulateModelFactoryFunction ffunc;
	try
	{
		// the next token is the name of the model
		// use "at" to check existence of the key in the models map
		ffunc = models.at(molstat::to_lower(tokens.front()));
	}
	catch(const out_of_range &e)
	{
		// model not found
		throw runtime_error("Model \"" + tokens.front() + "\" not found.");
	}
	molstat::SimulateModelFactory factory{ ffunc() };

	// start processing information from input...
	string command{ "not endmodel" };
	while(input && command != "endmodel")
	{
		string line;
		getline(input, line);

		// tokenize the string
		tokens = molstat::tokenize(line);
		if(tokens.size() == 0) // empty line
			continue;

		// the first token is the command name, pop it off and then process the
		// rest of the tokens
		command = molstat::to_lower(tokens.front());
		tokens.pop();

		// process the two available commands (get a submodel or set a
		// distribution)
		if(command == "distribution")
		{
			// get the name
			if(tokens.size() == 0) // no name specified
			{
				cout << "   No name specified on \"distribution\"." << endl;
			}
			else
			{
				const string name{ molstat::to_lower(tokens.front()) };
				tokens.pop();
				shared_ptr<const molstat::RandomDistribution> dist{ nullptr };

				// process the remaining tokens as a distribution
				try
				{
					dist = molstat::RandomDistributionFactory(move(tokens));

					// set the distribution
					factory.setDistribution(name, dist);
				}
				catch(const invalid_argument &e)
				{
					// print out a message, but don't error out
					cout << "   Error constructing random distribution: " <<
						e.what() << endl;
				}
			}
		}
		else if(command == "endmodel")
		{
			// don't need to do anything; we'll exit the input reading loop after
			// this cycle. we just don't want to print the error message.
		}
		else if(command == "model")
		{
			// construct a submodel... recursion!
		}
		else
		{
			cout << "   Unrecognized model command: \"" + command + "\"." << endl;
		}
	}

	// make sure the last command was endmodel; that is, we didn't hit EOF
	if(command != "endmodel")
		throw runtime_error("EOF encountered in \"model\" block.");

	return factory.getModel(); // let any exceptions pass up
}

molstat::ObservableIndex processObservable(
	std::queue<std::string> &&tokens,
	const std::map<std::string,
	          molstat::ObservableIndex> &observables)
{
	// make sure there is a token left, otherwise throw a runtime_error
	if(tokens.size() == 0)
		throw runtime_error("Observable name not specified.");

	molstat::ObservableIndex obsindex{ typeid(void*) };

	// look up the specified observable
	try
	{
		// the next token is the name of the observable
		// use "at" to check existence of the key in the observables map
		obsindex = observables.at(molstat::to_lower(tokens.front()));
	}
	catch(const out_of_range &e)
	{
		// observable not found
		throw runtime_error("Observable \"" + tokens.front() + "\" not found.");
	}

	return obsindex;
}
