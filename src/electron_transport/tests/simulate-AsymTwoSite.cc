/* This file is a part of MolStat, which is distributed under the Creative
   Commons Attribution-NonCommercial 4.0 International Public License.
   MolStat (c) 2014, Northwestern University. */
/**
 * \file tests/simulate-AsymTwoSite.cc
 * \brief Test suite for the asymmetric-coupling, two-site tight-binding
 *     model.
 *
 * \test Test suite for the asymmetric-coupling, two-site tight-binding
 *     model.
 *
 * \author Matthew G.\ Reuter
 * \date September 2014
 */

#include <cstdio>
#include <assert.h>
#include <cmath>
#include <vector>

#include "../simulator_models/asym_two_site_simulate_model.h"

using namespace std;

/**
 * \internal
 * \brief Main function for testing the asymmetric-coupling, two-site model.
 *
 * \param[in] argc The number of command-line arguments.
 * \param[in] argv The command-line arguments.
 * \return Exit status: 0 if the code passes the test, non-zero otherwise.
 * \endinternal
 */
int main(int argc, char **argv) {
	const double thresh = 1.0e-6;
	typedef molstat::AsymTwoSiteSimulateModel model;
	vector<double> params(6);

	// params[0] = EF, [1] = V, [2] = eps, [3] = gammal, [4] = gammar,
	//       [5] = beta

	// check known values for several parameter sets
	params = {0., 1., -4., 0.8, 1., -3.};
	assert(abs(0.121622 - model::transmission
		(params[0], 0., params[2], params[3], params[4], params[5])) < thresh);
	assert(abs(0.149936 - model::static_conductance(params)) < thresh);
	assert(abs(0.213248 - model::diff_conductance(params)) < thresh);

	params = {1., -0.4, -3., 0.4, 0.2, -0.8};
	assert(abs(0.000216257 - model::transmission
		(params[0], 0., params[2], params[3], params[4], params[5])) < thresh);
	assert(abs(0.000218231 - model::static_conductance(params)) < thresh);
	assert(abs(0.000222203 - model::diff_conductance(params)) < thresh);

	params = {-1., 1.4, 5., 0.67, 1.98, -1.6};
	assert(abs(0.00292927 - model::transmission
		(params[0], 0., params[2], params[3], params[4], params[5])) < thresh);
	assert(abs(0.00308371 - model::static_conductance(params)) < thresh);
	assert(abs(0.00340305 - model::diff_conductance(params)) < thresh);

	return 0;
}
