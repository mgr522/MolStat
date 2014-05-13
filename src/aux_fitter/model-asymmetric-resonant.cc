/*
This work is licensed under the Creative Commons Attribution 3.0 United States
License. To view a copy of this license, visit
http://creativecommons.org/licenses/by/3.0/us/ or send a letter to Creative
Commons, 444 Castro Street, Suite 900, Mountain View, California, 94041, USA.

Copyright (C) 2013 Oak Ridge National Laboratory
*/
/**
 * \file model-asymmetric-resonant.cc
 * \brief Implementation of the objective function and Jacobian for the
 *        asymmetric coupling, resonant tunneling model.
 *
 * This form is detailed in P.\ D.\ Williams, M.\ G.\ Reuter. J.\ Phys.\ 
 * Chem.\ C \b 117, 5937-5942 (2013).
 *
 * \author Matthew G.\ Reuter
 * \date August 2012, May 2013
 */

#include "models.h"
#include <gsl/gsl_math.h>

/// Struct for getting parameters into the integrals.
struct integral_data {
	double gamma1, gamma2, r, g;
};

// prototypes for the various model integrals
/**
 * \brief Integral for the asymmetric model objective function.
 *
 * \param[in] x The current value of x.
 * \param[in] params The parameters gamma1, gamma2, c, and g.
 * \return The value of the integrand.
 */
static double int_g(double x, void *params);

/**
 * \brief Integral for the asymmetric model objective function; part of the
 *        derivative with respect to r.
 *
 * \param[in] x The current value of x.
 * \param[in] params The parameters gamma1, gamma2, r, and g.
 * \return The value of the integrand.
 */
static double int_gr(double x, void *params);

/**
 * \brief Integral for the asymmetric model objective function; part of the
 *        derivative with respect to gamma1.
 *
 * \param[in] x The current value of x.
 * \param[in] params The parameters gamma1, gamma2, c, and g.
 * \return The value of the integrand.
 */
static double int_ggamma1(double x, void *params);

/**
 * \brief Integral for the asymmetric model objective function; part of the
 *        derivative with respect to gamma2.
 *
 * \param[in] x The current value of x.
 * \param[in] params The parameters gamma1, gamma2, c, and g.
 * \return The value of the integrand.
 */
static double int_ggamma2(double x, void *params);

int asymmetric_coupling_resonance_f(const gsl_vector *params, void *data,
	gsl_vector *resid) {

	integral_data id;
	double gi, pdfi, error, intmin, intmax;
	gsl_function f;
	f.function = &int_g;
	f.params = &id;

	// get the data values to be fitted
	size_t n = ((st_data*)data)->n;
	double *g = ((st_data*)data)->g;
	double *pdf = ((st_data*)data)->pdf;
	gsl_integration_workspace *w = ((st_data*)data)->w;

	// get the current parameters from GSL
	// this model uses gamma1, gamma2, and c
	id.gamma1 = gsl_vector_get(params, 0);
	id.gamma2 = gsl_vector_get(params, 1);
	id.r = gsl_vector_get(params, 2);
	double norm = gsl_vector_get(params, 3);
	
	size_t i;

	for(i = 0; i < n; ++i) {
		gi = g[i];
		id.g = gi;
		intmin = (2.0 - gi - 2.0*sqrt(1.0-gi)) / gi;
		intmax = (2.0 - gi + 2.0*sqrt(1.0-gi)) / gi;

		// calculate the integral
		gsl_integration_qags(&f, intmin, intmax, 0.0, 1.0e-7, 2000, w,
			&pdfi, &error);
		pdfi *= norm / (gi*sqrt(gi));
		gsl_vector_set(resid, i, pdfi - pdf[i]);
	}

	return GSL_SUCCESS;
}

int asymmetric_coupling_resonance_df(const gsl_vector *params,
	void *data, gsl_matrix *J) {

	integral_data id;
	double gi, dpdfdgamma1, dpdfdgamma2, dpdfdr, dpdfdnorm;
	double error, intmin, intmax;
	double intg, intgr, intggamma1, intggamma2;
	gsl_function f;
	f.params = &id;

	// get the data values to be fitted
	size_t n = ((st_data*)data)->n;
	double *g = ((st_data*)data)->g;
	gsl_integration_workspace *w = ((st_data*)data)->w;

	// get the current parameters from GSL
	// this model uses gamma1, gamma2, and c
	id.gamma1 = gsl_vector_get(params, 0);
	id.gamma2 = gsl_vector_get(params, 1);
	id.r = gsl_vector_get(params, 2);
	double norm = gsl_vector_get(params, 3);

	size_t i;

	for(i = 0; i < n; ++i) {
		gi = g[i];
		id.g = gi;
		intmin = (2.0 - gi - 2.0*sqrt(1.0-gi)) / gi;
		intmax = (2.0 - gi + 2.0*sqrt(1.0-gi)) / gi;

		// evaluate the four integrals
		f.function = &int_g;
		gsl_integration_qags(&f, intmin, intmax, 0.0, 1.0e-7, 2000, w,
			&intg, &error);

		f.function = &int_ggamma1;
		gsl_integration_qags(&f, intmin, intmax, 0.0, 1.0e-7, 2000, w,
			&intggamma1, &error);

		f.function = &int_ggamma2;
		gsl_integration_qags(&f, intmin, intmax, 0.0, 1.0e-7, 2000, w,
			&intggamma2, &error);

		f.function = &int_gr;
		gsl_integration_qags(&f, intmin, intmax, 0.0, 1.0e-7, 2000, w,
			&intgr, &error);

		// set the derivatives
		dpdfdgamma1 = norm / (gi*sqrt(gi)) * intggamma1;

		dpdfdgamma2 = norm / (gi*sqrt(gi)) * intggamma2;

		dpdfdr = -0.25 * norm * id.r * intgr * (id.gamma1*id.gamma1
			+ id.gamma2*id.gamma2) / (gi*gi*sqrt(gi));

		dpdfdnorm = norm * intg / (gi*sqrt(gi));

		gsl_matrix_set(J, i, 0, dpdfdgamma1);
		gsl_matrix_set(J, i, 1, dpdfdgamma2);
		gsl_matrix_set(J, i, 2, dpdfdr);
		gsl_matrix_set(J, i, 3, dpdfdnorm);
	}

	return GSL_SUCCESS;
}

int asymmetric_coupling_resonance_fdf(const gsl_vector *params, void *data,
	gsl_vector *resid, gsl_matrix *J) {

	// we can save some time since the int_g integral is needed for both
	// f and df...
	// combine the two functions above

	integral_data id;
	double gi, pdfi, dpdfdgamma1, dpdfdgamma2, dpdfdr, dpdfdnorm;
	double error, intmin, intmax;
	double intg, intgr, intggamma1, intggamma2;
	gsl_function f;
	f.params = &id;

	// get the data values to be fitted
	size_t n = ((st_data*)data)->n;
	double *g = ((st_data*)data)->g;
	double *pdf = ((st_data*)data)->pdf;
	gsl_integration_workspace *w = ((st_data*)data)->w;

	// get the current parameters from GSL
	// this model uses gamma1, gamma2, and c
	id.gamma1 = gsl_vector_get(params, 0);
	id.gamma2 = gsl_vector_get(params, 1);
	id.r = gsl_vector_get(params, 2);
	double norm = gsl_vector_get(params, 3);

	size_t i;

	for(i = 0; i < n; ++i) {
		gi = g[i];
		id.g = gi;
		intmin = (2.0 - gi - 2.0*sqrt(1.0-gi)) / gi;
		intmax = (2.0 - gi + 2.0*sqrt(1.0-gi)) / gi;

		// evaluate the four integrals
		f.function = &int_g;
		gsl_integration_qags(&f, intmin, intmax, 0.0, 1.0e-7, 2000, w,
			&intg, &error);

		f.function = &int_ggamma1;
		gsl_integration_qags(&f, intmin, intmax, 0.0, 1.0e-7, 2000, w,
			&intggamma1, &error);

		f.function = &int_ggamma2;
		gsl_integration_qags(&f, intmin, intmax, 0.0, 1.0e-7, 2000, w,
			&intggamma2, &error);

		f.function = &int_gr;
		gsl_integration_qags(&f, intmin, intmax, 0.0, 1.0e-7, 2000, w,
			&intgr, &error);

		// set the residual, can get dpdfdnorm cheaply in this process
		dpdfdnorm = intg / (gi*sqrt(gi));
		pdfi = norm * dpdfdnorm;

		// weight the data points, as described above
		gsl_vector_set(resid, i, pdfi - pdf[i]);

		// set the derivatives
		dpdfdgamma1 = norm / (gi*sqrt(gi)) * intggamma1;

		dpdfdgamma2 = norm / (gi*sqrt(gi)) * intggamma2;

		dpdfdr = -0.25 * id.r * intgr * (id.gamma1*id.gamma1
			+ id.gamma2*id.gamma2) / (gi*gi*sqrt(gi));

		gsl_matrix_set(J, i, 0, dpdfdgamma1);
		gsl_matrix_set(J, i, 1, dpdfdgamma2);
		gsl_matrix_set(J, i, 2, dpdfdr);
		gsl_matrix_set(J, i, 3, dpdfdnorm);
	}

	return GSL_SUCCESS;
}

static double int_g(double x, void *params) {
	integral_data &id = *(integral_data*)params;

	double temp1 = 4.0*x - id.g*(1.0+x)*(1.0+x);
	double temp2 = 1.0 + x*x;

	return x / (temp2 * sqrt(temp1 * temp2)) *
		(1.0 + (id.gamma1 + x*id.gamma2)*(id.gamma1 + x*id.gamma2) / temp2) *
		exp(-0.5*(x*id.gamma1 - id.gamma2)*(x*id.gamma1 - id.gamma2) / temp2
			-0.125*id.r*id.r*(id.gamma1*id.gamma1 + id.gamma2*id.gamma2)*temp1
			/ (temp2 * id.g));
}

static double int_gr(double x, void *params) {
	integral_data &id = *(integral_data*)params;

	double temp1 = 4.0*x - id.g*(1.0+x)*(1.0+x);
	double temp2 = 1.0 + x*x;

	return x*sqrt(temp1 / temp2) / (temp2 * temp2) *
		(1.0 + (id.gamma1 + x*id.gamma2)*(id.gamma1 + x*id.gamma2) / temp2) *
		exp(-0.5*(x*id.gamma1 - id.gamma2)*(x*id.gamma1 - id.gamma2) / temp2
			-0.125*id.r*id.r*(id.gamma1*id.gamma1 + id.gamma2*id.gamma2)*temp1
			/ (temp2 * id.g));
}

static double int_ggamma1(double x, void *params) {
	integral_data &id = *(integral_data*)params;

	double temp1 = 4.0*x - id.g*(1.0+x)*(1.0+x);
	double temp2 = 1.0 + x*x;
	double temp3 = id.gamma1 + x*id.gamma2;
	temp3 *= temp3 / temp2;

	return x / (temp2 * temp2 * sqrt(temp2)) *
		(((2.0-x*x)*id.gamma1 + 3.0*x*id.gamma2
			-temp3*x*(x*id.gamma1-id.gamma2)) / sqrt(temp1)
			-0.25*(1.0 + temp3)*id.r*id.r*id.gamma1*sqrt(temp1)/id.g) *
		exp(-0.5*(x*id.gamma1 - id.gamma2)*(x*id.gamma1 - id.gamma2) / temp2
			-0.125*id.r*id.r*(id.gamma1*id.gamma1 + id.gamma2*id.gamma2)*temp1
			/ (temp2 * id.g));
}

static double int_ggamma2(double x, void *params) {
	integral_data &id = *(integral_data*)params;

	double temp1 = 4.0*x - id.g*(1.0+x)*(1.0+x);
	double temp2 = 1.0 + x*x;
	double temp3 = id.gamma1 + x*id.gamma2;
	temp3 *= temp3 / temp2;

	return x / (temp2 * temp2 * sqrt(temp2)) *
		((3.0*x*id.gamma1 + (2.0*x*x-1.0)*id.gamma2
			+temp3*(x*id.gamma1-id.gamma2)) / sqrt(temp1)
			-0.25*(1.0 + temp3)*id.r*id.r*id.gamma2*sqrt(temp1)/id.g) *
		exp(-0.5*(x*id.gamma1 - id.gamma2)*(x*id.gamma1 - id.gamma2) / temp2
			-0.125*id.r*id.r*(id.gamma1*id.gamma1 + id.gamma2*id.gamma2)*temp1
			/ (temp2 * id.g));
}
