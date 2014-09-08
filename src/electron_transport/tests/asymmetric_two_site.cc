/* This file is a part of MolStat, which is distributed under the Creative
   Commons Attribution-NonCommercial 4.0 International Public License.
   MolStat (c) 2014, Northwestern University. */
/**
 * \file tests/asymmetric_two_site.cc
 * \brief Test suite for the asymmetric-coupling, two-site tight-binding
 *    model.
 *
 * \test Test suite for the asymmetric-coupling, two-site tight-binding
 *    model.
 *
 * \author Matthew G.\ Reuter
 * \date May 2014
 */

#include <cstdio>
#include <assert.h>
#include <cmath>
#include "../simulator_models/asymmetric_two_site.h"

using namespace std;

/**
 * \internal
 * \brief Main function for testing the asymmetric-coupling, one-site model.
 *
 * \param[in] argc The number of command-line arguments.
 * \param[in] argv The command-line arguments.
 * \return Exit status: 0 if the code passes the test, non-zero otherwise.
 * \endinternal
 */
int main(int argc, char **argv) {
	const double thresh = 1.0e-6;
	typedef AsymmetricTwoSiteModel model;

	// check known values for several parameter sets
	assert(abs(0.121622
		- model::transmission(0., -4., 0.8, 1.0, -3.)) < thresh);
	assert(abs(0.149936
		- model::static_conductance(0., 1., 0.5, -4., 0.8, 1.0, -3.)) < thresh);
	assert(abs(0.213248
		- model::diff_conductance(0., 1., 0.5, -4., 0.8, 1.0, -3.)) < thresh);

	assert(abs(0.000216257
		- model::transmission(1., -3., 0.4, 0.2, -0.8)) < thresh);
	assert(abs(0.000247897
		- model::static_conductance(1., -0.4, 0.8, -3., 0.4, 0.2, -0.8)) < thresh);
	assert(abs(0.000284847
		- model::diff_conductance(1., -0.4, 0.8, -3., 0.4, 0.2, -0.8)) < thresh);

	assert(abs(0.00292927
		- model::transmission(-1., 5., 0.67, 1.98, -1.6)) < thresh);
	assert(abs(0.00217835
		- model::static_conductance(-1., 1.4, 0.14, 5., 0.67, 1.98, -1.6)) < thresh);
	assert(abs(0.00164361
		- model::diff_conductance(-1., 1.4, 0.14, 5., 0.67, 1.98, -1.6)) < thresh);

	return 0;
}
