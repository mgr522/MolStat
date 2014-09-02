/* This file is a part of MolStat, which is distributed under the Creative
   Commons Attribution-NonCommercial 4.0 International Public License.
   MolStat (c) 2014, Northwestern University. */
/**
 * \file tests/asymmetric_voltage_one_site.cc
 * \brief Test suite for the voltage-dependent, asymmetric-coupling, single-site
 *    tight-binding model.
 *
 * \test Test suite for the asymmetric-coupling, single-site tight-binding
 *    model for conductances.
 *
 * \author Matthew G.\ Reuter
 * \date September 2014
 */

#include <cstdio>
#include <assert.h>
#include <cmath>

#include "../simulator_models/asymmetric_voltage_one_site.h"

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
	typedef AsymmetricVoltageOneSiteModel model;

	// check known values for several parameter sets
	assert(abs(0.0499376
		- model::transmission(0., 1., -4., 0.8, 1.0, 0.1)) < thresh);
	assert(abs(0.0506743
		- model::static_conductance(0., 1., 0.5, -4., 0.8, 1.0, 0.1)) < thresh);
	assert(abs(0.0546687
		- model::diff_conductance(0., 1., 0.5, -4., 0.8, 1.0, 0.1)) < thresh);

	assert(abs(0.000801201
		- model::transmission(1., -0.4, -9., 0.4, 0.2, -0.03)) < thresh);
	assert(abs(0.000821124
		- model::static_conductance(1., -0.4, 0.8, -9., 0.4, 0.2, -0.03)) <
		thresh);
	assert(abs(0.000843753
		- model::diff_conductance(1., -0.4, 0.8, -9., 0.4, 0.2, -0.03)) <
		thresh);

	assert(abs(0.00399256
		- model::transmission(3., 1.4, -17., 0.67, 1.98, 1.3)) < thresh);
	assert(abs(0.00422874
		- model::static_conductance(3., 1.4, 0.14, -17., 0.67, 1.98, 1.3)) <
		thresh);
	assert(abs(0.00534931
		- model::diff_conductance(3., 1.4, 0.14, -17., 0.67, 1.98, 1.3)) <
		thresh);

	return 0;
}
