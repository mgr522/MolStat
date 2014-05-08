/**
 * \file symmetric_one_site.cc
 * \brief Test suite for the symmetric-coupling, single-site tight-binding
 *     model.
 *
 * \author Matthew G.\ Reuter
 * \date May 2014
 */

#include <cstdio>
#include <assert.h>
#include <cmath>
#include "../aux_simulator/symmetric_one_site.h"

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
	typedef SymmetricOneSiteModel model;

	// check known values for several parameter sets
	assert(abs(0.0384615
		- model::transmission(0., -4., 0.8)) < thresh);
	assert(abs(0.0390172
		- model::static_conductance(0., 1., 0.5, -4., 0.8)) < thresh);
	assert(abs(0.0401438
		- model::diff_conductance(0., 1., 0.5, -4., 0.8)) < thresh);

	assert(abs(0.00159744
		- model::transmission(1., -9., 0.4)) < thresh);
	assert(abs(0.00163709
		- model::static_conductance(1., -0.4, 0.8, -9., 0.4)) < thresh);
	assert(abs(0.00167814
		- model::diff_conductance(1., -0.4, 0.8, -9., 0.4)) < thresh);

	assert(abs(0.00112099
		- model::transmission(3., -17., 0.67)) < thresh);
	assert(abs(0.00118115
		- model::static_conductance(3., 1.4, 0.14, -17., 0.67)) < thresh);
	assert(abs(0.00124526
		- model::diff_conductance(3., 1.4, 0.14, -17., 0.67)) < thresh);

	return 0;
}
