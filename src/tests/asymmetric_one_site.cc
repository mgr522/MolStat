/**
 * \file tests/asymmetric_one_site.cc
 * \brief Test suite for the asymmetric-coupling, single-site tight-binding
 *     model.
 *
 * \author Matthew G.\ Reuter
 * \date May 2014
 */

#include <cstdio>
#include <assert.h>
#include <cmath>
#include "../aux_simulator/asymmetric_one_site.h"

using namespace std;

/**
 * \brief Main function for testing the symmetric-coupling, one-site model.
 *
 * \param[in] argc The number of command-line arguments.
 * \param[in] argv The command-line arguments.
 * \return Exit status: 0 if the code passes the test, non-zero otherwise.
 */
int main(int argc, char **argv) {
	const double thresh = 1.0e-6;
	typedef AsymmetricOneSiteModel model;

	// check known values for several parameter sets
	assert(abs(0.0475907
		- model::transmission(0., -4., 0.8, 1.0)) < thresh);
	assert(abs(0.0482617
		- model::static_conductance(0., 1., 0.5, -4., 0.8, 1.0)) < thresh);
	assert(abs(0.0496212
		- model::diff_conductance(0., 1., 0.5, -4., 0.8, 1.0)) < thresh);

	assert(abs(0.000799281
		- model::transmission(1., -9., 0.4, 0.2)) < thresh);
	assert(abs(0.000819131
		- model::static_conductance(1., -0.4, 0.8, -9., 0.4, 0.2)) < thresh);
	assert(abs(0.000839689
		- model::diff_conductance(1., -0.4, 0.8, -9., 0.4, 0.2)) < thresh);

	assert(abs(0.00330201
		- model::transmission(3., -17., 0.67, 1.98)) < thresh);
	assert(abs(0.00347858
		- model::static_conductance(3., 1.4, 0.14, -17., 0.67, 1.98)) < thresh);
	assert(abs(0.00366672
		- model::diff_conductance(3., 1.4, 0.14, -17., 0.67, 1.98)) < thresh);

	return 0;
}
