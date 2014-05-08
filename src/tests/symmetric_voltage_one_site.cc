/**
 * \file tests/symmetric_voltage_one_site.cc
 * \brief Test suite for the symmetric-coupling, voltage-dependent, single-site
 *     tight-binding model.
 *
 * \author Matthew G.\ Reuter
 * \date May 2014
 */

#include <cstdio>
#include <assert.h>
#include <cmath>
#include "../aux_simulator/symmetric_voltage_one_site.h"

using namespace std;

/**
 * \brief Main function for testing the symmetric-coupling, voltage-
 *    dependent, one-site model.
 *
 * \param[in] argc The number of command-line arguments.
 * \param[in] argv The command-line arguments.
 * \return Exit status: 0 if the code passes the test, non-zero otherwise.
 */
int main(int argc, char **argv) {
	const double thresh = 1.0e-6;
	typedef SymmetricVoltageOneSiteModel model;

	// check known values for several parameter sets
	assert(abs(0.06639
		- model::transmission(0., 1., -4., 0.8)) < thresh);
	assert(abs(0.0679934
		- model::static_conductance(0., 1., 0.5, -4., 0.8)) < thresh);
	assert(abs(0.114507
		- model::diff_conductance(0., 1., 0.5, -4., 0.8)) < thresh);

	assert(abs(0.00147710
		- model::transmission(1., -0.4, -9., 0.4)) < thresh);
	assert(abs(0.00151231
		- model::static_conductance(1., -0.4, 0.8, -9., 0.4)) < thresh);
	assert(abs(0.00143116
		- model::diff_conductance(1., -0.4, 0.8, -9., 0.4)) < thresh);

	assert(abs(0.00129587
		- model::transmission(3., 1.4, -17., 0.67)) < thresh);
	assert(abs(0.00137100
		- model::static_conductance(3., 1.4, 0.14, -17., 0.67)) < thresh);
	assert(abs(0.00166364
		- model::diff_conductance(3., 1.4, 0.14, -17., 0.67)) < thresh);

	return 0;
}
