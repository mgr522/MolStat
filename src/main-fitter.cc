/* This file is a part of MolStat, which is distributed under the Creative
   Commons Attribution-NonCommercial 4.0 International Public License.
   MolStat (c) 2014, Northwestern University. */
/**
 * \file main-fitter.cc
 * \brief Main function for fitting single-molecule data to the desired
 *    functional form.
 *
 * This program fits data (using one of the implemented models) and outputs the
 * best-fit parameters. Because all of these functional forms are non-linear,
 * multiple initial guesses are used, and the overall best fit is output.
 *
 * \author Matthew G.\ Reuter
 * \date September 2014
 */

#include <memory>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <list>
#include <utility>
#include <map>
#include <string>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_multifit_nlin.h>

#include "general/string_tools.h"
#include "general/histogram_tools/bin_style.h"
#include "electron_transport/fitter_models/transport_fit_models.h"

using namespace std;

/**
 * \internal
 * \brief Main function.
 *
 * Parses the input parameters and outputs the best-fit parameters (and,
 * optionally, iteration-to-iteration details).
 *
 * \param[in] argc The number of command-line arguments.
 * \param[in] argv The command-line arguments.
 * \return Exit status; 0 for normal.
 * \endinternal
 */
int main(int argc, char **argv) {
	// non-linear fitting tools
	unique_ptr<gsl_multifit_fdfsolver,
	           decltype(&gsl_multifit_fdfsolver_free)>
		solver(nullptr, &gsl_multifit_fdfsolver_free);
	gsl_multifit_function_fdf fdf;

	// status of the fit
	unique_ptr<gsl_vector, decltype(&gsl_vector_free)>
		vec(nullptr, &gsl_vector_free);
	bool hasfit, iterprint, usedefaultguess;
	double resid, bestresid;
	int status;

	// the model we're fitting against
	unique_ptr<molstat::FitModel<1>> model;
	// the list of models:
	// stored as a map of string (model name) to a model factory, given the
	// list of data
	map<string, molstat::FitModelFactory<1>> models;

	// the data we're fitting
	list<pair<array<double, 1>, double>> data;

	// various parameters for determining the best fit
	vector<double> bestfit;
	list<vector<double>> initvals;
	list<vector<double>>::const_iterator initval;
	unique_ptr<molstat::BinStyle> binstyle;

	// counters & auxiliary variables
	size_t i, iter;
	double g, pdf;
	FILE *f;
	string line, modelname;
	vector<string> tokens;

	gsl_set_error_handler_off();

	// load the models
	// FitModelAdd calls appear here
	molstat::load_transport_models(models);

	// set up the fit -- read in parameters from stdin
	// Line 1: One token specifying the model to use for fitting
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
	// instantiating the model)
	modelname = tokens[0];

	// Line 2: The file name of the conductance histogram data to fit
	try {
		line = molstat::getline(stdin);
	} catch(const runtime_error &e) {
		fprintf(stderr, "Error: %s.\n", e.what());
		return 0;
	}
	molstat::tokenize(line, tokens);
	if(tokens.size() < 1) {
		fprintf(stderr, "Error: file name expected in line 2.\n");
		return 0;
	}

	// read in the data points from the specified file
	f = fopen(tokens[0].c_str(), "r");
	if(!f) {
		fprintf(stderr, "Error opening %s for input.\n", tokens[0].c_str());
		return 0;
	}
	// read all lines in the file.
	while(fscanf(f, "%le %le", &g, &pdf) != EOF) {
		data.emplace_back(pair<array<double, 1>, double>({{g}}, pdf));
	}
	fclose(f);
	size_t nbin = data.size();

	// set the model
	try {
		// models.at() returns a function for instantiating the model
		model = models.at(molstat::to_lower(modelname))(data);
	}
	catch(const out_of_range &e) {
		fprintf(stderr, "Error: model \"%s\" not found.\n", modelname.c_str());
		return 0;
	}
	// set the vector size for the number of fitting parameters
	bestfit.resize(model->nfit);

	// set up GSL details
	fdf = model->gsl_handle();
	vec.reset(gsl_vector_alloc(model->nfit));

	// Remaining lines: auxiliary options
	// default options
	iterprint = false; // don't print details at every iteration
	initvals.clear(); // no initial guesses
	usedefaultguess = false; // user specifies to use the default guesses
	binstyle = molstat::BinStyleFactory({{"linear"}}); // default: linear bins

	// process the lines and override any of the default options
	try {
		while(true) {
			line = molstat::getline(stdin);

			molstat::tokenize(line, tokens);
			if(tokens.size() > 0) {
				line = tokens[0];
				molstat::make_lower(line);

				if(line == "print")
					iterprint = true;
				else if(line == "noprint")
					iterprint = false;
				else if(line == "guess") {
					// this line specifies a guess -- is it to use the defaults
					// or a specific set of parameters?
					if(tokens.size() == 1) {
						fprintf(stderr, "Error: No initial guess specified. " \
							"Skipping line.\n");
					}
					else {
						line = tokens[1];
						molstat::make_lower(line);
						if(line == "default")
							usedefaultguess = true;
						else {
							// this is a user-specified initial guess

							// remove the first token from the vector ("guess")
							tokens.erase(tokens.begin());

							// process the initial guess
							try {
								model->append_initial_guess(tokens, initvals);
							}
							catch(const invalid_argument &e) {
								fprintf(stderr, "Error: %s Skipping input line.\n",
									e.what());
							}
						}
					}
				}
				else if(line == "bin"){
					// load a bin style

					// first need to remove the "bin" token
					tokens.erase(tokens.begin());
					molstat::make_lower(tokens[0]);

					try {
						binstyle = molstat::BinStyleFactory(tokens);
					}
					catch(const invalid_argument &e) {
						fprintf(stderr, "Error: Unknown binning style. Skipping " \
							"line.\n");
					}
				}
				// add other keywords/options here
			}
		}
	}
	catch(const runtime_error &e) {
		// this just means we hit EOF -- stop trying to read more
	}

	// use the bin type to "unmask", if necessary, the data so that we fit in
	// g, not some function of g
	for(list<pair<array<double, 1>, double>>::iterator iter = data.begin();
		iter != data.end(); ++iter) {

		// transform the independent variable back to g
		(*iter).first[0] = binstyle->invmask((*iter).first[0]);

		// transform the PDF: P_g(g)
		//    = P_mask(g)(mask(g)) * dmaskdx(invmask(u))
		(*iter).second *= binstyle->dmaskdx((*iter).first[0]);
	}

	// do we need to load the default guesses?
	if(usedefaultguess || initvals.size() == 0)
		model->append_default_guesses(initvals);

	solver.reset(
		gsl_multifit_fdfsolver_alloc(gsl_multifit_fdfsolver_lmsder, nbin,
			model->nfit));

	// we don't have a successful fit at the start
	hasfit = false;

	// perform fits with all the initial values
	for(initval = initvals.cbegin(); initval != initvals.cend(); ++initval) {
		// load the initial values
		for(i = 0; i < model->nfit; ++i)
			gsl_vector_set(vec.get(), i, (*initval)[i]);

		gsl_multifit_fdfsolver_set(solver.get(), &fdf, vec.get());

		// start iterating
		iter = 0;
		if(iterprint) {
			fprintf(stdout, "Iter=%3zu, ", iter);
			model->print_fit(stdout, molstat::gsl_to_std(solver->x));
			fprintf(stdout, "\n");
		}

		do {
			++iter;
			status = gsl_multifit_fdfsolver_iterate(solver.get());
			if(iterprint) {
				fprintf(stdout, "Iter=%3zu, ", iter);
				model->print_fit(stdout, molstat::gsl_to_std(solver->x));
				fprintf(stdout, "\n");
			}

			if(status)
				break;

			status = gsl_multifit_test_delta(solver->dx, solver->x, 1.0e-4,
				1.0e-4);
		} while(status == GSL_CONTINUE && iter < 1000);

		if(iterprint) {
			if(status != GSL_CONTINUE && status != GSL_SUCCESS)
				printf("   %s\n\n", gsl_strerror(status));
		}

		// did we converge, iteration out, or error out?
		if(status != GSL_CONTINUE && status != GSL_SUCCESS &&
			status != GSL_ENOPROG) {

			// errored out
			continue; // try the next guess of initial values
		}

		// do some final processing
		resid = gsl_blas_dnrm2(solver->f);
		if(iterprint)
			printf("Residual = %.6e\n\n", resid);

		if(!hasfit || resid < bestresid) {
			// copy over the parameters
			bestresid = resid;
			for(i = 0; i < model->nfit; ++i)
				bestfit[i] = gsl_vector_get(solver->x, i);

			hasfit = true;
		}
	}

	// did we get a fit?
	if(!hasfit) {
		fprintf(stderr, "Error fitting.\n");
	}
	else {
		// make sure the fit parameters are good
		model->process_fit_parameters(bestfit);

		// print ou the fit
		printf("Resid = %.6e\n", bestresid);
		model->print_fit(stdout, bestfit);
		printf("\n");
	}

	return 0;
}
