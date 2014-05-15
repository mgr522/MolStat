/*
This work is licensed under the Creative Commons Attribution 3.0 United States
License. To view a copy of this license, visit
http://creativecommons.org/licenses/by/3.0/us/ or send a letter to Creative
Commons, 444 Castro Street, Suite 900, Mountain View, California, 94041, USA.

Copyright (C) 2013 Oak Ridge National Laboratory
*/
/**
 * \file model-symmetric-resonant.cc
 * \brief Implementation of the objective function and Jacobian for the
 *        symmetric-coupling, resonant tunneling model.
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

int symmetric_coupling_resonance_f(const gsl_vector *params, void *data,
	gsl_vector *resid) {

	double gi, pdfi;

	// get the data values to be fitted
	size_t n = ((st_data*)data)->n;
	double *g = ((st_data*)data)->g;
	double *pdf = ((st_data*)data)->pdf;

	// get the current parameters from GSL
	// this model uses gamma
	double gamma = gsl_vector_get(params, 0);
	double norm = gsl_vector_get(params, 1);
	
	size_t i;

	for(i = 0; i < n; ++i) {
		gi = g[i];
		pdfi = norm / sqrt(gi*gi*gi*(1.0 - gi))
			* exp(-0.5*gamma*gamma*(1.0 - gi) / gi);

		// owing to the singularity in the form -- the data can span several
		// orders of magnitude with most points much smaller than a few --
		// scale by the size of the point to give more weight to the smaller
		// points
		gsl_vector_set(resid, i, (pdfi - pdf[i]) / pdf[i]);
	}

	return GSL_SUCCESS;
}

int symmetric_coupling_resonance_df(const gsl_vector *params,
	void *data, gsl_matrix *J) {

	double gi, dpdfdgamma, dpdfdnorm;

	// get the data values to be fitted
	size_t n = ((st_data*)data)->n;
	double *g = ((st_data*)data)->g;
	double *pdf = ((st_data*)data)->pdf;

	// get the current parameters from GSL
	// this model uses c
	double gamma = gsl_vector_get(params, 0);
	double norm = gsl_vector_get(params, 1);

	size_t i;

	for(i = 0; i < n; ++i) {
		gi = g[i];
		dpdfdgamma = -gamma * norm * sqrt((1.0 - gi) / gi)
			* exp(-0.5*gamma*gamma*(1.0-gi)/gi) / (gi*gi);
		dpdfdgamma /= pdf[i]; // scaling as described above

		dpdfdnorm = 1.0 / sqrt(gi*gi*gi*(1.0 - gi))
			* exp(-0.5*gamma*gamma*(1.0 - gi) / gi);
		dpdfdnorm /= pdf[i]; // scaling, again

		gsl_matrix_set(J, i, 0, dpdfdgamma);
		gsl_matrix_set(J, i, 1, dpdfdnorm);
	}

	return GSL_SUCCESS;
}

int symmetric_coupling_resonance_fdf(const gsl_vector *params, void *data,
	gsl_vector *resid, gsl_matrix *J) {

	symmetric_coupling_resonance_f(params, data, resid);
	symmetric_coupling_resonance_df(params, data, J);

	return GSL_SUCCESS;
}
