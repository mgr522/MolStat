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
 * \date May 2014
 */

#include <gsl/gsl_blas.h>

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <list>
#include <utility>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_multifit_nlin.h>
#include "aux_fitter/fit_model_1d.h"

using namespace std;

#if 0
/**
 * \brief Prints the status of the nonlinear fit after each iteration.
 *
 * \param[in] iter The present iteration.
 * \param[in] s The GSL solver.
 * \param[in] model The model being used.
 */
void iteration_print(size_t iter, shared_ptr<gsl_multifit_fdfsolver> s,
	char model);
#endif

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
	// non-linear fitting tools
	shared_ptr<gsl_multifit_fdfsolver> solver;
	gsl_multifit_function_fdf fdf;

	// status of the fit
	shared_ptr<gsl_vector> vec, bestfit, beststdev;
	bool hasfit, iterprint;
	double resid, bestresid;
	int status;

	// model & parameters
	shared_ptr<FitModel<1>> model;
	list<pair<array<double, 1>, double>> data;
	list<vector<double>> initvals;
	list<vector<double>>::const_iterator initval;

	// counters & auxiliary variables
	size_t i, iter;

	gsl_set_error_handler_off();

	// read in the data points from the STDIN
	size_t nbin = 39;
	for(i = 0; i < nbin; ++i) {
		double g, pdf;
		scanf("%le %le", &g, &pdf);

		data.emplace_back(pair<array<double, 1>, double>({{g}}, pdf));
	}

	// set the model
	model = get_fit_model("symmetricresonant", data);

	// set up GSL details
	fdf = model->gsl_handle();
	vec.reset(gsl_vector_alloc(model->nfit), &gsl_vector_free);
	bestfit.reset(gsl_vector_alloc(model->nfit), &gsl_vector_free);
	beststdev.reset(gsl_vector_alloc(model->nfit), &gsl_vector_free);

#if 0
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

		// allocate numerical integration workspace
		data.w = gsl_integration_workspace_alloc(2000);
	}
	normparams = params + 1; // number of parameters, including the norm
#endif

	solver.reset(
		gsl_multifit_fdfsolver_alloc(gsl_multifit_fdfsolver_lmsder, nbin,
			model->nfit),
		&gsl_multifit_fdfsolver_free);

	// we don't have a successful fit at the start
	hasfit = false;

	// perform fits with all the initial values
	initvals = model->initial_guesses();
	for(initval = initvals.cbegin(); initval != initvals.cend(); ++initval) {
		// load the initial values
		for(i = 0; i < model->nfit; ++i) {
			gsl_vector_set(vec.get(), i, (*initval)[i]);
		}
		gsl_multifit_fdfsolver_set(solver.get(), &fdf, vec.get());

		// start iterating
		iter = 0;
#if 0
		if(iterprint)
			iteration_print(iter, solver, model);
#endif

		do {
			++iter;
			status = gsl_multifit_fdfsolver_iterate(solver.get());
#if 0
			if(iterprint)
				iteration_print(iter, solver, model);
#endif

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
			// copy over the parameters
			bestresid = resid;
			for(i = 0; i < model->nfit; ++i)
				gsl_vector_set(bestfit.get(), i, gsl_vector_get(solver->x, i));

			hasfit = true;
		}
	}

	// did we get a fit?
	if(!hasfit) {
		fprintf(stderr, "Error fitting.\n");
	}
#if 0
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
#endif

	return 0;
}

#if 0
void iteration_print(size_t iter, shared_ptr<gsl_multifit_fdfsolver> s,
	char model) {

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
#endif
