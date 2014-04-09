/*
This work is licensed under the Creative Commons Attribution 3.0 United States
License. To view a copy of this license, visit
http://creativecommons.org/licenses/by/3.0/us/ or send a letter to Creative
Commons, 444 Castro Street, Suite 900, Mountain View, California, 94041, USA.

Copyright (C) 2013 Oak Ridge National Laboratory
*/
/**
 * \file model-symmetric-nonresonant.cc
 * \brief Implementation of the objective function and Jacobian for the
 *        symmetric-coupling, nonresonant tunneling model.
 *
 * This form is detailed in P.\ D.\ Williams, M.\ G.\ Reuter. J.\ Phys.\ 
 * Chem.\ C \b 117, 5937-5942 (2013).
 *
 * \author Patrick D.\ Williams and Matthew G.\ Reuter
 * \date July 2012, May 2013
 */

#include "models.h"
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_erf.h>

int symmetric_coupling_nonresonance_f(const gsl_vector *params, void *data,
	gsl_vector *resid) {

	double gi, pdfi;

	// get the data values to be fitted
	size_t n = ((st_data*)data)->n;
	double *g = ((st_data*)data)->g;
	double *pdf = ((st_data*)data)->pdf;

	// get the current parameters from GSL
	// this model uses c and d
	double c = gsl_vector_get(params, 0);
	double d = gsl_vector_get(params, 1);
	double norm = gsl_vector_get(params, 2);

	size_t i;

	for(i = 0; i < n; ++i) {
		gi = g[i];
		double cd1 = c*sqrt(gi) - d*sqrt(1.0 - gi);
		pdfi = norm / ((1.0-gi)*sqrt(gi*(1.0-gi)))
			* exp(-0.5*cd1*cd1/(1.0-gi));
		gsl_vector_set(resid, i, pdfi - pdf[i]);
	}

	return GSL_SUCCESS;
}

int symmetric_coupling_nonresonance_df(const gsl_vector *params,
	void *data, gsl_matrix *J) {

	double gi, dpdfdc, dpdfdd, dpdfdnorm;

	// get the data values to be fitted
	size_t n = ((st_data*)data)->n;
	double *g = ((st_data*)data)->g;

	// get the current parameters from GSL
	// this model uses c and d
	double c = gsl_vector_get(params, 0);
	double d = gsl_vector_get(params, 1);
	double norm = gsl_vector_get(params, 2);

	size_t i;

	for(i = 0; i < n; ++i) {
		gi = g[i];
		double cd1 = c*sqrt(gi) - d*sqrt(1.0 - gi);
		double expcd = exp(-0.5*cd1*cd1 / (1.0 - gi));

		dpdfdc = -norm * cd1 * expcd / ((1.0-gi)*(1.0-gi)*sqrt(1.0-gi));

		dpdfdd = norm * cd1 * expcd / ((1.0-gi)*(1.0-gi)*sqrt(gi));

		dpdfdnorm = expcd / ((1.0-gi)*sqrt(gi*(1.0-gi)));

		gsl_matrix_set(J, i, 0, dpdfdc);
		gsl_matrix_set(J, i, 1, dpdfdd);
		gsl_matrix_set(J, i, 2, dpdfdnorm);
	}

	return GSL_SUCCESS;
}

int symmetric_coupling_nonresonance_fdf(const gsl_vector *params,
	void *data, gsl_vector *resid, gsl_matrix *J) {

	symmetric_coupling_nonresonance_f(params, data, resid);
	symmetric_coupling_nonresonance_df(params, data, J);

	return GSL_SUCCESS;
}
