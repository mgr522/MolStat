/**
 * \file test_gsl_std_vector.cc
 * \brief Test suite for interconverting gsl and std vector types.
 *
 * \author Matthew G.\ Reuter
 * \date May 2014
 */

#include <cstdio>
#include <assert.h>
#include "fit_model_interface.h"

using namespace std;

/**
 * \brief Main function for testing the Histogram1D class with linear binning.
 *
 * \param[in] argc The number of command-line arguments.
 * \param[in] argv The command-line arguments.
 * \return Exit status: 0 if the code passes the test, non-zero otherwise.
 */
int main(int argc, char **argv) {
	size_t i;
	const size_t size = 10;
	gsl_vector *gslv;
	vector<double> vec;

	gslv = gsl_vector_alloc(size);

	for(i = 0; i < size; ++i)
		gsl_vector_set(gslv, i, 1.*i);

	vec = gsl_to_std(gslv);

	for(i = 0; i < size; ++i)
		assert(abs(vec[i] - gsl_vector_get(gslv, i)) < 1.0e-6);

	return 0;
}
