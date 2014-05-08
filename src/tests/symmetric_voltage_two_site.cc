/**
 * \file tests/symmetric_voltage_two_site.cc
 * \brief Test suite for the symmetric-coupling, voltage-dependent, two-site
 *     tight-binding model.
 *
 * \author Matthew G.\ Reuter
 * \date May 2014
 */

#include <cstdio>
#include <assert.h>
#include <cmath>
#include "../aux_simulator/symmetric_voltage_two_site.h"

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
	typedef SymmetricVoltageTwoSiteModel model;

	// check known values for several parameter sets
	assert(abs(0.107326
		- model::transmission(0., 1., -4., 0.8, -3.)) < thresh);
	assert(abs(0.136415
		- model::static_conductance(0., 1., 0.5, -4., 0.8, -3.)) < thresh);
	assert(abs(0.222768
		- model::diff_conductance(0., 1., 0.5, -4., 0.8, -3.)) < thresh);

	assert(abs(0.000433828
		- model::transmission(1., -0.4, -3., 0.4, -0.8)) < thresh);
	assert(abs(0.000497406
		- model::static_conductance(1., -0.4, 0.8, -3., 0.4, -0.8)) < thresh);
	assert(abs(0.000577222
		- model::diff_conductance(1., -0.4, 0.8, -3., 0.4, -0.8)) < thresh);

	assert(abs(0.631062
		- model::transmission(3., 1.4, 1.1, 0.67, -1.6)) < thresh);
	assert(abs(0.457519
		- model::static_conductance(3., 1.4, 0.14, 1.1, 0.67, -1.6)) < thresh);
	assert(abs(-0.0108761
		- model::diff_conductance(3., 1.4, 0.14, 1.1, 0.67, -1.6)) < thresh);

	return 0;
}
