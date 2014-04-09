/*
This work is licensed under the Creative Commons Attribution 3.0 United States
License. To view a copy of this license, visit
http://creativecommons.org/licenses/by/3.0/us/ or send a letter to Creative
Commons, 444 Castro Street, Suite 900, Mountain View, California, 94041, USA.

Copyright (C) 2013 Oak Ridge National Laboratory
*/
/**
 * \file main-fitter.cc
 * \brief Main function for fitting histogram data to the desired functional
 *        form.
 *
 * This program fits conductance data (using one of the models listed below)
 * and outputs the best-fit parameters. Since all of these functional forms
 * are non-linear, multiple initial guesses are used, and the overall best fit
 * is output.
 *
 * Available functional forms are:
 *    - `n' The non-resonant tunneling (symmetric channel-lead coupling) model.
 *    - `r' The resonant tunneling (symmetric channel-lead coupling) model.
 *    - `a' The resonant tunneling (asymmetric channel-lead coupling) model.
 *
 * There are two required command-line arguments, and another optional one.
 *    -# The number of bins in the histogram (number of data points to read
 *       from standard input).
 *    -# The model to use (see above for character flags).
 *    -# [Optional] Produce iteration-to-iteration output for the nonlinear
 *       fitting process. This can produce a lot of output.
 *
 * \author Matthew G.\ Reuter
 * \date July 2012, May 2013
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_multifit_nlin.h>
#include <queue>

#include "models.h"

using namespace std;

/**
 * \brief Prints the status of the nonlinear fit after each iteration.
 *
 * \param[in] iter The present iteration.
 * \param[in] solver The GSL solver.
 * \param[in] model The model being used.
 */
void iteration_print(size_t iter, gsl_multifit_fdfsolver *s, char model);

/**
 * \brief Main function.
 *
 * Parses the input parameters and outputs the best-fit parameters (and,
 * optionally, iteration-to-iteration details).
 *
 * \param[in] argc The number of command-line arguments.
 * \param[in] argv The command-line arguments.
 * \return Exit status; 0 for normal.
 */
int main(int argc, char **argv) {
	int status;
	size_t i, iter, nbin;
	size_t params, normparams;
	double *g, *pdf;
	gsl_multifit_fdfsolver *solver;
	gsl_multifit_function_fdf fdf;
	st_data data;
	gsl_vector *vec, *bestfit, *beststdev;
	gsl_matrix *covar;
	double chi, dof, rel;
	char model;
	queue<double> initvals;
	bool hasfit, iterprint;
	double resid, bestresid;
	double fits[3];
	double norm;

	// get the command-line arguments
	if(argc != 3 && argc != 4) {
		fprintf(stderr, "Usage: ./fitter nbin model [print]\n" \
			"   nbin is the number of bins in the histogram.\n" \
			"   model is the model to fit against.\n" \
			"      - 'n' for non-resonant tunneling, symmetric coupling\n" \
			"      - 'r' for resonant tunneling, symmetric coupling\n" \
			"      - 'a' for resonant tunneling, asymmetric coupling\n" \
			"   [print] is an optional argument that produces output after each" \
			" iteration of the\n           non-linear fit routine. This can "\
			"produce a LOT of output.\n" \
			"NOTE: The data is expected through stdin.\n");
		return 0;
	}

	nbin = atoi(argv[1]);
	if(nbin < 1) {
		fprintf(stderr, "Error: Use at least one bin.\n");
		return 0;
	}

	model = argv[2][0];
	if(model != 'r' && model != 'n' && model != 'a') {
		fprintf(stderr, "Error: Model must be 'r', 'n', or 'a' for the\n" \
			"resonant (symmetric), non-resonant (symmetric), and resonant\n" \
			"(asymmetric) models, respectively.\n");
		return 0;
	}

	iterprint = (argc == 4);

	gsl_set_error_handler_off();

	g = (double*)malloc(nbin*sizeof(double));
	pdf = (double*)malloc(nbin*sizeof(double));

	// read in the data points from STDIN
	for(i = 0; i < nbin; ++i)
		scanf("%le %le", g + i, pdf + i);

	// setup the data for GSL
	data.n = nbin;
	data.g = g;
	data.pdf = pdf;
	data.w = NULL;
	fdf.n = nbin;
	fdf.params = &data;

	// a note on the initial values:
	// since each model uses a various number of parameters, the initvals
	// queue just stores each parameter as an item.
	// For example, on each trial in the 'r' model, only 1 item is read
	// from initvals; but for each trial in 'n', 2 items are read.

	// model-specific setup
	if(model == 'n') {
		params = 2; // number of parameters in the fit

		// set up a range for initial values
		double cvals[] = {50.0, 100.0, 200.0, 300.0, 400.0, 500.0};
		double dvals[] = {5.0, 10.0, 20.0, 30.0, 40.0, 50.0};
		for(int i = 0; i < 6; ++i)
		for(int j = 0; j < 6; ++j) {
			initvals.push(cvals[i]);
			initvals.push(dvals[j]);
		}

		// initialize GSL
		fdf.f = &symmetric_coupling_nonresonance_f;
		fdf.df = &symmetric_coupling_nonresonance_df;
		fdf.fdf = &symmetric_coupling_nonresonance_fdf;
	}
	else if(model == 'r') {
		params = 1; // number of parameters in the fit

		// set up some initial value guesses
		initvals.push(5.0);
		initvals.push(10.0);
		initvals.push(20.0);
		initvals.push(35.0);
		initvals.push(50.0);

		// initialize GSL
		fdf.f = &symmetric_coupling_resonance_f;
		fdf.df = &symmetric_coupling_resonance_df;
		fdf.fdf = &symmetric_coupling_resonance_fdf;
	}
	else if(model == 'a') {
		params = 3; // number of parameters in the fit

		// set up a range for initial values
		double gammavals[] = {5.0, 10.0, 20.0, 30.0, 40.0};
		double cvals[] = {0.1, 0.5, 1.0, 2.0, 10.0};
		for(int i = 0; i < 5; ++i)
		for(int j = 0; j <= i; ++j)
		for(int k = 0; k < 5; ++k) {
			initvals.push(gammavals[i]);
			initvals.push(gammavals[j]);
			initvals.push(cvals[k]);
		}

		// initialize GSL
		fdf.f = &asymmetric_coupling_resonance_f;
		fdf.df = &asymmetric_coupling_resonance_df;
		fdf.fdf = &asymmetric_coupling_resonance_fdf;

		// allocate numerical integration workspace
		data.w = gsl_integration_workspace_alloc(2000);
	}
	normparams = params + 1; // number of parameters, including the norm

	vec = gsl_vector_alloc(normparams);
	bestfit = gsl_vector_alloc(normparams);
	beststdev = gsl_vector_alloc(normparams);
	fdf.p = normparams;

	solver = gsl_multifit_fdfsolver_alloc(gsl_multifit_fdfsolver_lmsder, nbin,
		normparams);

	// we don't have a successful fit at the start
	hasfit = false;

	// go through all the initial values
	do {
		// load the initial values
		for(i = 0; i < params; ++i) {
			gsl_vector_set(vec, i, initvals.front());
			initvals.pop();
		}
		gsl_vector_set(vec, params, 1.0);
		gsl_multifit_fdfsolver_set(solver, &fdf, vec);

		// start iterating
		iter = 0;
		if(iterprint)
			iteration_print(iter, solver, model);

		do {
			++iter;
			status = gsl_multifit_fdfsolver_iterate(solver);
			if(iterprint)
				iteration_print(iter, solver, model);

			if(status)
				break;

			status = gsl_multifit_test_delta(solver->dx, solver->x, 1.0e-4,
				1.0e-4);
		} while(status == GSL_CONTINUE && iter < 1000);

		if(iterprint) {
			if(status != GSL_CONTINUE && status != GSL_SUCCESS) {
				printf("   %s\n\n", gsl_strerror(status));
			}
			else {
				printf("\n");
			}
		}

		// did we converge, iteration out, or error out?
		if(status != GSL_CONTINUE && status != GSL_SUCCESS &&
			status != GSL_ENOPROG) {

			// errored out
			continue; // try the next guess of initial values
		}

		// do some final processing
		resid = gsl_blas_dnrm2(solver->f);
		if(!hasfit || resid < bestresid) {
			covar = gsl_matrix_alloc(normparams, normparams);
			gsl_multifit_covar(solver->J, 0.0, covar);
			chi = gsl_blas_dnrm2(solver->f);
			dof = nbin - normparams;
			rel = GSL_MAX_DBL(1, chi / sqrt(dof));

			// copy over the parameters and standard deviations
			bestresid = resid;
			for(i = 0; i < normparams; ++i) {
				gsl_vector_set(bestfit, i, gsl_vector_get(solver->x, i));
				gsl_vector_set(beststdev, i,
					rel*sqrt(gsl_matrix_get(covar, i, i)));
			}

			gsl_matrix_free(covar);
			hasfit = true;
		}
	} while(!initvals.empty());

	// did we get a fit?
	if(!hasfit) {
		fprintf(stderr, "Error fitting.\n");
	}
	else {
		// do some post-processing (i.e., the fit might have let the gammas
		// become negative, etc.)
		if(model == 'n') {
			// no problems seen, as yet
		}
		else if(model == 'r') {
			fits[0] = gsl_vector_get(bestfit, 0);
			norm = gsl_vector_get(bestfit, 1);

			// gamma might be negative, fix it
			if(fits[0] < 0.0) {
				gsl_vector_set(bestfit, 0, -fits[0]);
			}
		}
		else if(model == 'a') {
			fits[0] = gsl_vector_get(bestfit, 0);
			fits[1] = gsl_vector_get(bestfit, 1);
			fits[2] = gsl_vector_get(bestfit, 2);
			norm = gsl_vector_get(bestfit, 3);

			// both gammas might be negative, fix it
			if(fits[0] < 0.0 && fits[1] < 0.0) {
				gsl_vector_set(bestfit, 0, -fits[0]);
				gsl_vector_set(bestfit, 1, -fits[1]);
			}

			// r might be negative, fix it
			if(fits[2] < 0.0) {
				gsl_vector_set(bestfit, 2, -fits[2]);
			}
		}

		// print ou the fit
		printf("Resid = %.6e\n", bestresid);
		printf("Constant of proportionality = %.6e\n",
			gsl_vector_get(bestfit, params));
		if(model == 'n') {
			printf("c = %.6e\n", gsl_vector_get(bestfit, 0));
			printf("d = %.6e\n", gsl_vector_get(bestfit, 1));
		}
		else if(model == 'r') {
			printf("gamma = %.6e\n", gsl_vector_get(bestfit, 0));
		}
		else if(model == 'a') {
			printf("gammaL = %.6e\n", gsl_vector_get(bestfit, 0));
			printf("gammaR = %.6e\n", gsl_vector_get(bestfit, 1));
			printf("r = %.6e\n", gsl_vector_get(bestfit, 2));
		}
	}

	// free all allocations
	if(model == 'a')
		gsl_integration_workspace_free(data.w);
	gsl_multifit_fdfsolver_free(solver);
	gsl_vector_free(beststdev);
	gsl_vector_free(bestfit);
	gsl_vector_free(vec);
	free(pdf);
	free(g);

	return 0;
}

void iteration_print(size_t iter, gsl_multifit_fdfsolver *s, char model) {
	if(model == 'n') {
		printf("iter %04u: c=%.4e d=%.4e resid=%.4e\n", (unsigned int)iter,
			gsl_vector_get(s->x, 0), gsl_vector_get(s->x, 1),
			gsl_blas_dnrm2(s->f));
	}
	else if(model == 'r') {
		printf("iter %04u: c=%.4e resid=%.4e\n", (unsigned int)iter,
			gsl_vector_get(s->x, 0), gsl_blas_dnrm2(s->f));
	}
	else if(model == 'a') {
		printf("iter %04u: gamma1=%.4e gamma2=%.4e c=%.4e resid=%.4e\n",
			(unsigned int)iter, gsl_vector_get(s->x, 0), gsl_vector_get(s->x, 1),
			gsl_vector_get(s->x, 2), gsl_blas_dnrm2(s->f));
	}
}
