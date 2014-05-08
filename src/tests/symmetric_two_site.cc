/**
 * \file tests/symmetric_two_site.cc
 * \brief Test suite for the symmetric-coupling, two-site
 *     tight-binding model.
 *
 * \author Matthew G.\ Reuter
 * \date May 2014
 */

#include <cstdio>
#include <assert.h>
#include <cmath>
#include "../aux_simulator/symmetric_two_site.h"

using namespace std;

/**
 * \brief Main function for testing the symmetric-coupling, two-site model.
 *
 * \param[in] argc The number of command-line arguments.
 * \param[in] argv The command-line arguments.
 * \return Exit status: 0 if the code passes the test, non-zero otherwise.
 */
int main(int argc, char **argv) {
	const double thresh = 1.0e-6;
	typedef SymmetricTwoSiteModel model;

	// check known values for several parameter sets
	assert(abs(0.101007
		- model::transmission(0., 1., -4., 0.8, -3.)) < thresh);
	assert(abs(0.127042
		- model::static_conductance(0., 1., 0.5, -4., 0.8, -3.)) < thresh);
	assert(abs(0.186815
		- model::diff_conductance(0., 1., 0.5, -4., 0.8, -3.)) < thresh);

	assert(abs(0.000431590
		- model::transmission(1., -0.4, -3., 0.4, -0.8)) < thresh);
	assert(abs(0.000494645
		- model::static_conductance(1., -0.4, 0.8, -3., 0.4, -0.8)) < thresh);
	assert(abs(0.000568266
		- model::diff_conductance(1., -0.4, 0.8, -3., 0.4, -0.8)) < thresh);

	assert(abs(0.459683
		- model::transmission(3., 1.4, 1.1, 0.67, -1.6)) < thresh);
	assert(abs(0.560554
		- model::static_conductance(3., 1.4, 0.14, 1.1, 0.67, -1.6)) < thresh);
	assert(abs(0.230111
		- model::diff_conductance(3., 1.4, 0.14, 1.1, 0.67, -1.6)) < thresh);

	return 0;
}
