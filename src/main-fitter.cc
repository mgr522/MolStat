/* This file is a part of MolStat, which is distributed under the Creative
   Commons Attribution-NonCommercial 4.0 International Public License.

   (c) 2014 Northwestern University. */

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
#include <fstream>
#include <iostream>
#include <iomanip>
#include <stdexcept>
#include <cstdlib>
#include <cmath>
#include <list>
#include <utility>
#include <map>
#include <string>
#include <limits>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_multifit_nlin.h>

#include <config.h>

#include "general/string_tools.h"
#include "general/histogram_tools/bin_style.h"
#include "general/histogram_tools/bin_linear.h"
#include "general/fitter_tools/fit_model_interface.h"

#if BUILD_TRANSPORT_FITTER
#include "electron_transport/fitter_models/transport_fit_module.h"
#endif

using namespace std;

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
int main(int argc, char **argv)
{
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

	// counters & auxiliary variables
	size_t i, iter, maxiter;
	string line, modelname;

	gsl_set_error_handler_off();

	// load the models
	// FitModelAdd calls appear here
	#if BUILD_TRANSPORT_FITTER
	molstat::transport::load_models(models);
	#endif

	// set up the fit -- read in parameters from stdin
	// Line 1: One token specifying the model to use for fitting
	if(cin)
		getline(cin, line);
	else
	{
		cerr << "Error: EOF encountered in line 1.\n" << endl;
		return 0;
	}

	molstat::TokenContainer tokens = molstat::tokenize(line);
	if(tokens.size() < 1)
	{
		cerr << "Error: model name expected in line 1." << endl;
		return 0;
	}
	// store the model name for later (need to continue reading data before
	// instantiating the model)
	modelname = tokens.front();

	// Line 2: The file name of the conductance histogram data to fit
	if(cin)
		getline(cin, line);
	else
	{
		cerr << "Error: EOF encountered in line 2." << endl;
		return 0;
	}
	tokens = molstat::tokenize(line);
	if(tokens.size() < 1)
	{
		cerr << "Error: file name expected in line 2." << endl;
		return 0;
	}

	// read in the data points from the specified file
	{
		ifstream f(tokens.front());
		if(!f)
		{
			cerr << "Error opening " << tokens.front() << " for input." << endl;
			return 0;
		}
		// read all lines in the file.
		while(f)
		{
			double g, pdf;
			f >> g >> pdf;
			data.emplace_back(pair<array<double, 1>, double>({{g}}, pdf));
		}
		f.close();
	}
	size_t nbin = data.size();

	// set the model
	try
	{
		// models.at() returns a function for instantiating the model
		model = models.at(molstat::to_lower(modelname))(data);
	}
	catch(const out_of_range &e)
	{
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
	maxiter = 100; // only allow 100 iterations per initial guess
	initvals.clear(); // no initial guesses
	usedefaultguess = false; // user specifies to use the default guesses
	// default is linear bins
	unique_ptr<molstat::BinStyle> binstyle{ new molstat::BinLinear(1) };

	// process the lines and override any of the default options
	try
	{
		while(getline(cin, line)) // read until EOF
		{
			tokens = molstat::tokenize(line);
			if(tokens.size() > 0)
			{
				line = molstat::to_lower(tokens.front());
				tokens.pop();

				if(line == "print")
					iterprint = true;
				else if(line == "noprint")
					iterprint = false;
				else if(line == "guess")
				{
					// this line specifies a guess -- is it to use the defaults
					// or a specific set of parameters?
					if(tokens.size() == 0)
					{
						cerr << "Error: No initial guess specified. Skipping line."
							<< endl;
					}
					else
					{
						line = molstat::to_lower(tokens.front());

						if(line == "default")
						{
							tokens.pop();
							usedefaultguess = true;
						}
						else // this is a user-specified initial guess
						{
							// process the initial guess
							try
							{
								model->append_initial_guess(move(tokens), initvals);
							}
							catch(const invalid_argument &e)
							{
								cerr << "Error: " << e.what() <<
									" Skipping input line." << endl;
							}
						}
					}
				}
				else if(line == "bin") // load a bin style
				{
					if(tokens.size() == 0)
					{
						cerr << "Error: No binning style specified. Skipping line."
							<< endl;
					}
					else
					{
						// need to create an auxiliary TokenContainer with the number
						// of bins (not needed) in the front spot
						molstat::TokenContainer tc;
						tc.push("1");
						tc.push(molstat::to_lower(tokens.front()));
						tokens.pop();
						while(!tokens.empty())
						{
							tc.push(tokens.front());
							tokens.pop();
						}

						try
						{
							binstyle = molstat::BinStyleFactory(move(tc));
						}
						catch(const invalid_argument &e)
						{
							cerr << "Error: Unknown binning style. Skipping line." <<
								endl;
						}
					}
				}
				else if(line == "maxiter") // change the max number of iterations
				{
					if(tokens.size() == 0)
					{
						cerr << "Error: Number of max iterations unspecified." \
							" Skipping line." << endl;
					}
					else
					{
						try
						{
							maxiter = molstat::cast_string<size_t>(tokens.front());
							tokens.pop();
						}
						catch(const bad_cast &e)
						{
							cerr << "Error interpreting number of max iterations." \
								" Skipping line." << endl;
						}
					}
				}
				// add other keywords/options here
			}
		}
	}
	catch(const runtime_error &e)
	{
		// this just means we hit EOF -- stop trying to read more
	}

	// use the bin type to "unmask", if necessary, the data so that we fit in
	// g, not some function of g
	for(list<pair<array<double, 1>, double>>::iterator iter = data.begin();
		iter != data.end(); ++iter)
	{

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

	// we don't have a successful fit at the start... set the residual as high
	// as possible
	hasfit = false;
	bestresid = std::numeric_limits<double>::max();

	// perform fits with all the initial values
	for(initval = initvals.cbegin(); initval != initvals.cend(); ++initval)
	{
		// load the initial values
		for(i = 0; i < model->nfit; ++i)
			gsl_vector_set(vec.get(), i, (*initval)[i]);

		gsl_multifit_fdfsolver_set(solver.get(), &fdf, vec.get());

		// start iterating
		iter = 0;
		if(iterprint)
		{
			cout << "Iter=" << setw(3) << iter << ", ";
			model->print_fit(cout, molstat::gsl_to_std(solver->x));
			cout << endl;
		}

		do
		{
			++iter;
			status = gsl_multifit_fdfsolver_iterate(solver.get());
			if(iterprint)
			{
				cout << "Iter=" << setw(3) << iter << ", ";
				model->print_fit(cout, molstat::gsl_to_std(solver->x));
				cout << endl;
			}

			if(status)
				break;

			status = gsl_multifit_test_delta(solver->dx, solver->x, 1.0e-4,
				1.0e-4);
		} while(status == GSL_CONTINUE && iter < maxiter);

		if(iterprint)
		{
			if(status != GSL_CONTINUE && status != GSL_SUCCESS)
				cout << "   " << gsl_strerror(status) << '\n' << endl;
		}

		// did we converge, iteration out, or error out?
		if(status != GSL_CONTINUE && status != GSL_SUCCESS &&
			status != GSL_ENOPROG)
		{
			// errored out
			continue; // try the next guess of initial values
		}

		// do some final processing
		resid = gsl_blas_dnrm2(solver->f);
		if(iterprint)
			cout << "Residual = " << scientific << setprecision(6) << resid <<
				'\n' << endl;

		if(!hasfit || resid < bestresid)
		{
			// copy over the parameters
			bestresid = resid;
			for(i = 0; i < model->nfit; ++i)
				bestfit[i] = gsl_vector_get(solver->x, i);

			hasfit = true;
		}
	}

	// did we get a fit?
	if(!hasfit)
	{
		cerr << "Error fitting." << endl;
	}
	else
	{
		// make sure the fit parameters are good
		model->process_fit_parameters(bestfit);

		// print ou the fit
		cout << "Resid = " << scientific << setprecision(6) << bestresid << '\n';
		model->print_fit(cout, bestfit);
		cout << endl;
	}

	return 0;
}
