/* This file is a part of MolStat, which is distributed under the Creative
   Commons Attribution-NonCommercial 4.0 International Public License.
   MolStat (c) 2014, Northwestern University. */
/**
 * \file tests/simulate-SymTwoSite.cc
 * \brief Test suite for the symmetric-coupling, two-site tight-binding
 *     model.
 *
 * \test Test suite for the symmetric-coupling, two-site tight-binding
 *     model.
 *
 * \author Matthew G.\ Reuter
 * \date September 2014
 */

#include <cstdio>
#include <assert.h>
#include <cmath>
#include <vector>

#include "../simulator_models/sym_two_site_simulate_model.h"

using namespace std;

/**
 * \internal
 * \brief Main function for testing the symmetric-coupling, two-site model.
 *
 * \param[in] argc The number of command-line arguments.
 * \param[in] argv The command-line arguments.
 * \return Exit status: 0 if the code passes the test, non-zero otherwise.
 * \endinternal
 */
int main(int argc, char **argv) {
	const double thresh = 1.0e-6;
	typedef SymTwoSiteSimulateModel model;
	vector<double> params(5);

	// params[0] = EF, [1] = V, [2] = eps, [3] = gamma, [4] = beta

	// check known values for several parameter sets
	params = {0., 1., -4., 0.8, -3.};
	assert(abs(0.101007 - model::transmission
		(params[0], 0., params[2], params[3], params[4])) < thresh);
	assert(abs(0.127042 - model::static_conductance(params)) < thresh);
	assert(abs(0.186815 - model::diff_conductance(params)) < thresh);

	params = {1., -0.4, -3., 0.4, -0.8};
	assert(abs(0.000431590 - model::transmission
		(params[0], 0., params[2], params[3], params[4])) < thresh);
	assert(abs(0.000435520 - model::static_conductance(params)) < thresh);
	assert(abs(0.000443426 - model::diff_conductance(params)) < thresh);

	params = {3., 1.4, 1.1, 0.67, -1.6};
	assert(abs(0.459683 - model::transmission
		(params[0], 0., params[2], params[3], params[4])) < thresh);
	assert(abs(0.480791 - model::static_conductance(params)) < thresh);
	assert(abs(0.294527 - model::diff_conductance(params)) < thresh);

	return 0;
}
