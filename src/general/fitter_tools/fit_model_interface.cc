/**
 * \file fit_model_interface.cc
 * \brief Implementations for non-templated functions in fit_model.h.
 *
 * \author Matthew G.\ Reuter
 * \date May 2014
 */

#include "fit_model_interface.h"

std::vector<double> gsl_to_std(const gsl_vector *gslv) {
	std::vector<double> ret(gslv->size);

	for(std::size_t i = 0; i < gslv->size; ++i)
		ret[i] = gsl_vector_get(gslv, i);

	return ret;
}
