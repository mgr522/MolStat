/*
This work is licensed under the Creative Commons Attribution 3.0 United States
License. To view a copy of this license, visit
http://creativecommons.org/licenses/by/3.0/us/ or send a letter to Creative
Commons, 444 Castro Street, Suite 900, Mountain View, California, 94041, USA.

Copyright (C) 2013 Oak Ridge National Laboratory
*/
/**
 * \file models.h
 * \brief Defines data structures for interfacing the conductance histogram
 *        models with the GSL non-linear fitting methods.
 *
 * Also includes each model's header for ease later.
 *
 * \author Matthew G.\ Reuter
 * \date May 2013
 */

#ifndef __models_h__
#define __models_h__

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_integration.h>

/**
 * \brief Struct for passing conductance and histogram data to the various
 *        residual/Jacobian GSL functions.
 */
typedef struct {
	/// The number of data points to be fit.
	size_t n;

	/// The conductance values.
	double *g;

	/// The histogram values [estimated PDF(g)].
	double *pdf;

	/// The workspace for numerical integration.
	gsl_integration_workspace *w;
} st_data;

#include "model-symmetric-nonresonant.h"
#include "model-symmetric-resonant.h"
#include "model-asymmetric-resonant.h"

#endif
