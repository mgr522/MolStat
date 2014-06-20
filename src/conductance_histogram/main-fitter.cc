/* This file is a part of MolStat, which is distributed under the Creative
   Commons Attribution-NonCommercial 4.0 International Public License.
   MolStat (c) 2014, Northwestern University. */
/**
 * \file main-fitter.cc
 * \brief Main function for fitting histogram data to the desired functional
 *        form.
 *
 * This program fits conductance data (using one of the implemented models)
 * and outputs the best-fit parameters. Since all of these functional forms
 * are non-linear, multiple initial guesses are used, and the overall best fit
 * is output.
 *
 * \author Matthew G.\ Reuter
 * \date June 2014
 */

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <list>
#include <utility>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_multifit_nlin.h>

#include <general/string_tools.h>
#include "fitter_models/cond_hist_fit_model.h"

using namespace std;

/**
 * \brief Main function.
 *
 * Parses the input parameters and outputs the best-fit parameters (and,
 * optionally, iteration-to-iteration details).
 *
 * \todo Add options for bin types.
 *
 * \param[in] argc The number of command-line arguments.
 * \param[in] argv The command-line arguments.
 * \return Exit status; 0 for normal.
 */
int main(int argc, char **argv) {
	// non-linear fitting tools
	shared_ptr<gsl_multifit_fdfsolver> solver;
	gsl_multifit_function_fdf fdf;

	// status of the fit
	shared_ptr<gsl_vector> vec;
	bool hasfit, iterprint, usedefaultguess;
	double resid, bestresid;
	int status;

	// model & parameters
	shared_ptr<FitModel<1>> model;
	list<pair<array<double, 1>, double>> data;
	vector<double> bestfit;
	list<vector<double>> initvals;
	list<vector<double>>::const_iterator initval;

	// counters & auxiliary variables
	size_t i, iter;
	double g, pdf;
	FILE *f;
	string line, modelname;
	vector<string> tokens;

	gsl_set_error_handler_off();

	// set up the fit -- read in parameters from stdin
	// Line 1: One token specifying the model to use for fitting
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
	// invoking the constructor
	modelname = tokens[0];
	make_lower(modelname);

	// Line 2: The file name of the conductance histogram data to fit
	try {
		line = getline(stdin);
	} catch(const runtime_error &e) {
		fprintf(stderr, "Error: %s.\n", e.what());
		return 0;
	}
	tokenize(line, tokens);
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
		model = get_cond_hist_fit_model(modelname, data);
	} catch(const invalid_argument &e) {
		fprintf(stderr, "Error: unknown model '%s'.\n", modelname.c_str());
		return 0;
	}
	bestfit.resize(model->nfit);

	// set up GSL details
	fdf = model->gsl_handle();
	vec.reset(gsl_vector_alloc(model->nfit), &gsl_vector_free);

	// Remaining lines: auxiliary options
	// default options
	iterprint = false; // don't print details at every iteration
	initvals.clear(); // no initial guesses
	usedefaultguess = false; // user specifies to use the default guesses

	// process the lines and override any of the default options
	try {
		while(true) {
			line = getline(stdin);

			tokenize(line, tokens);
			if(tokens.size() > 0) {
				line = tokens[0];
				make_lower(line);

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
						make_lower(line);
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
								fprintf(stderr, "ERROR: %s\n", e.what());
							}
						}
					}
				}
				// add other keywords/options here
			}
		}
	}
	catch(const runtime_error &e) {
		// this just means we hit EOF -- stop trying to read more
	}

	// do we need to load the default guesses?
	if(usedefaultguess || initvals.size() == 0)
		model->append_default_guesses(initvals);

	solver.reset(
		gsl_multifit_fdfsolver_alloc(gsl_multifit_fdfsolver_lmsder, nbin,
			model->nfit),
		&gsl_multifit_fdfsolver_free);

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
			model->print_fit(stdout, gsl_to_std(solver->x));
			fprintf(stdout, "\n");
		}

		do {
			++iter;
			status = gsl_multifit_fdfsolver_iterate(solver.get());
			if(iterprint) {
				fprintf(stdout, "Iter=%3zu, ", iter);
				model->print_fit(stdout, gsl_to_std(solver->x));
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
