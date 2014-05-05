/**
 * \file main-simulator.cc
 * \brief Main function for simulating conductance data using Landauer theory.
 *
 * This code reads in various input parameters from standard in, and uses these
 * parameters to simulate conductance data. The conductance data is designed to
 * be binned into conductance histograms. Both zero-bias (1D) and voltage-
 * dependent (2D) conductance data/histograms can be simulated.
 *
 * \author Matthew G.\ Reuter
 * \date May 2014
 */

#include <memory>
#include <cstdio>
#include <stdexcept>
#include <string>
#include <vector>
#include <functional>
#include <map>
#include "string_tools.h"
#include "aux_random_distributions/rng.h"
#include "aux_simulator/symmetric_one_site.h"
#include "aux_simulator/asymmetric_one_site.h"
#include "aux_simulator/symmetric_voltage_one_site.h"

using namespace std;

/**
 * \brief Enum for the types of conductance calculations we have.
 */
enum class CalculationType {
	Static, // Voltage-Dependent Static Conductance
	Differential, // Voltage-Dependend Differential Conductance
	ZeroBias // Zero-Bias Differential Conductance
};

/**
 * \brief Takes the list of distributions from the input deck and constructs
 *    a model of the specified type.
 *
 * \throw invalid_argument if any required distributions are not specified.
 *
 * \param[in] str The name of the type, assumed to be lowercase.
 * \param[in] parameters Map of distribution names (lowercase) to the
 *    distributions.
 * \return A shared pointer to the model object.
 */
shared_ptr<ConductanceModel> set_distributions(const string str,
	const map<string, shared_ptr<RandomDistribution>> &parameters);

/**
 * \brief Main function for simulating a histogram.
 *
 * Parses the input parameters and outputs randomly-generated conductance
 * data.
 *
 * \param[in] argc The number of command-line arguments.
 * \param[in] argv The command-line arguments.
 * \return Exit status; 0 for normal.
 */
int main(int argc, char **argv) {

	// initialize the GSL random number generator
	gsl_rng_env_setup();
	shared_ptr<gsl_rng> r(gsl_rng_alloc(gsl_rng_default), &gsl_rng_free);
	gsl_rng_set(r.get(), 0xFEEDFACE);

	CalculationType type;
	int i, n;
	double EF;
	map<string, shared_ptr<RandomDistribution>> parameters;
	shared_ptr<ConductanceModel> model;
	shared_ptr<RandomDistribution> dist_V, dist_eta;
	function<void(void)> conductance_function;

	string line, modeltype, name;
	vector<string> tokens;

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

	// make a link to the correct model creator (need to continue reading data,
	// store for later)
	modeltype = tokens[0];
	make_lower(modeltype);

	// Line 2: Type of conductance to calculate
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
	if(tokens[0] == "static")
		type = CalculationType::Static;
	else if(tokens[0] == "differential")
		type = CalculationType::Differential;
	else if(tokens[0] == "zerobias")
		type = CalculationType::ZeroBias;
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
		n = stoi(tokens[0]);
	}
	catch(const invalid_argument &e) {
		fprintf(stderr, "Error: unrecognizable number '%s'.\n",
			tokens[0].c_str());
		return 0;
	}

	// Line 4: Fermi level
	try {
		line = getline(stdin);
	}
	catch(const runtime_error &e) {
		fprintf(stderr, "Error: %s.\n", e.what());
		return 0;
	}
	tokenize(line, tokens);
	if(tokens.size() < 1) {
		fprintf(stderr, "Error Fermi energy expected in line 4.\n");
		return 0;
	}
	try {
		EF = stod(tokens[0]);
	}
	catch(const invalid_argument &e) {
		fprintf(stderr, "Error: unable to parse %s to the Fermi energy.\n",
			tokens[0].c_str());
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
					parameters[name] = distribution_from_tokens(tokens);
				}
				catch(const invalid_argument &e) {
					fprintf(stderr, "Error: unable to form a random number " \
						"distribution from:\n   %s\n%s\n", line.c_str(), e.what());
				}
			}
		}
	}
	catch(const runtime_error &e) {
	}

	// set the distributions required by the model
	try {
		model = set_distributions(modeltype, parameters);
	}
	catch(const invalid_argument &e) {
		fprintf(stderr, "Error initializing the model: %s\n", e.what());
		return 0;
	}

	// set up the run function for each type...
	// this may require other distributions to be set
	if(type == CalculationType::Static ||
		type == CalculationType::Differential) {

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

		// set the conductance function
		if(type == CalculationType::Static) {
			conductance_function = [=] () {
				double V = dist_V->sample(r);
				double GV = model->static_conductance(r, EF, dist_eta->sample(r),
					V);
				printf("%.6f %.6f\n", V, GV);
			};
		}
		else if(type == CalculationType::Differential) {
			conductance_function = [=] () {
				double V = dist_V->sample(r);
				double GV = model->diff_conductance(r, EF, dist_eta->sample(r), V);
				printf("%.6f %.6f\n", V, GV);
			};
		}
	}
	else if(type == CalculationType::ZeroBias) {
		// no extra distributions required... proceed to the calculation
		conductance_function = [=] () {
			printf("%.6f\n", model->zero_bias_conductance(r, EF));
		};
	}
	else {
		// should never be here
		fprintf(stderr, "An unknown error has occured.\n");
		return 0;
	}

	// Get the requested number of voltage-conductance data points
	for (i = 0; i < n; ++i) {
		conductance_function();
	}

	return 0;
}

shared_ptr<ConductanceModel> set_distributions(const string str,
	const map<string, shared_ptr<RandomDistribution>> &parameters) {

	// create the actual model and process the random number distributions
	if(str == "symmetriconesitemodel") {
		shared_ptr<RandomDistribution> dist_gamma, dist_eps;

		// populate the gamma and epsilon distributions
		try {
			dist_gamma = parameters.at("gamma");
		}
		catch(const out_of_range &e) {
			throw invalid_argument("A distribution for \"gamma\" must be " \
				"specified.");
		}

		try {
			dist_eps = parameters.at("epsilon");
		}
		catch(const out_of_range &e) {
			throw invalid_argument("A distribution for \"epsilon\" must be " \
				"specified.");
		}

		return make_shared<SymmetricOneSiteModel>(dist_eps, dist_gamma);
	}
	else if(str == "symmetricvoltageonesitemodel") {
		shared_ptr<RandomDistribution> dist_gamma, dist_eps;

		// populate the gamma and epsilon distributions
		try {
			dist_gamma = parameters.at("gamma");
		}
		catch(const out_of_range &e) {
			throw invalid_argument("A distribution for \"gamma\" must be " \
				"specified.");
		}

		try {
			dist_eps = parameters.at("epsilon");
		}
		catch(const out_of_range &e) {
			throw invalid_argument("A distribution for \"epsilon\" must be " \
				"specified.");
		}

		return make_shared<SymmetricVoltageOneSiteModel>(dist_eps, dist_gamma);
	}
	else if(str == "asymmetriconesitemodel") {
		shared_ptr<RandomDistribution> dist_gammaL, dist_gammaR, dist_eps;

		// populate the gamma and epsilon distributions
		try {
			dist_gammaL = parameters.at("gammal");
		}
		catch(const out_of_range &e) {
			throw invalid_argument("A distribution for \"gammaL\" must be " \
				"specified.");
		}

		try {
			dist_gammaR = parameters.at("gammar");
		}
		catch(const out_of_range &e) {
			throw invalid_argument("A distribution for \"gammaR\" must be " \
				"specified.");
		}

		try {
			dist_eps = parameters.at("epsilon");
		}
		catch(const out_of_range &e) {
			throw invalid_argument("A distribution for \"epsilon\" must be " \
				"specified.");
		}

		return make_shared<AsymmetricOneSiteModel>(dist_eps, dist_gammaL,
			dist_gammaR);
	}
	else
		throw invalid_argument("Unrecognized model. Options are:\n" \
			"   SymmetricOneSiteModel - " \
				"Symmetric-Coupling, One-Site Model\n" \
			"   SymmetricVoltageOneSiteModel - " \
				"Symmetric-Coupling, Voltage-Dependent One-Site Model\n" \
			"   AsymmetricOneSiteModel - " \
				"Asymmetric-Coupling, One-Site Model\n");

	// should never be here
	throw invalid_argument("Shouldn't be here.");
	return nullptr;
}
