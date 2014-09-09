/* This file is a part of MolStat, which is distributed under the Creative
   Commons Attribution-NonCommercial 4.0 International Public License.
   MolStat (c) 2014, Northwestern University. */
/**
 * \file tests/simulate-AsymOneSite.cc
 * \brief Test suite for the asymmetric-coupling, single-site tight-binding
 *     model.
 *
 * \test Test suite for the asymmetric-coupling, single-site tight-binding
 *     model.
 *
 * \author Matthew G.\ Reuter
 * \date September 2014
 */

#include <cstdio>
#include <assert.h>
#include <cmath>
#include <vector>

#include "../simulator_models/asym_one_site_simulate_model.h"

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
	typedef AsymOneSiteSimulateModel model;
	vector<double> params(6);

	// params[0] = EF, [1] = V, [2] = eps, [3] = gammal, [4] = gammar, [5] = a

	// check known values for several parameter sets
	params = {0., 1., -4., 0.8, 1., 0.};
	assert(abs(0.0475907 - model::transmission
		(params[0], 0., params[2], params[3], params[4], params[5])) < thresh);
	assert(abs(0.0482617 - model::static_conductance(params)) < thresh);
	assert(abs(0.0496212 - model::diff_conductance(params)) < thresh);

	params = {1., -0.4, -9., 0.4, 0.2, 0.};
	assert(abs(0.000799281 - model::transmission
		(params[0], 0., params[2], params[3], params[4], params[5])) < thresh);
	assert(abs(0.000799600 - model::static_conductance(params)) < thresh);
	assert(abs(0.000800238 - model::diff_conductance(params)) < thresh);

	params = {3., 1.4, -17., 0.67, 1.98, 0.};
	assert(abs(0.00330201 - model::transmission
		(params[0], 0., params[2], params[3], params[4], params[5])) < thresh);
	assert(abs(0.00330602 - model::static_conductance(params)) < thresh);
	assert(abs(0.00331404 - model::diff_conductance(params)) < thresh);

	params = {0., 1., -4., 0.8, 1., 0.1};
	assert(abs(0.0475907 - model::transmission
		(params[0], 0., params[2], params[3], params[4], params[5])) < thresh);
	assert(abs(0.0506743 - model::static_conductance(params)) < thresh);
	assert(abs(0.0546687 - model::diff_conductance(params)) < thresh);

	params = {1., -0.4, -9., 0.4, 0.2, -0.03};
	assert(abs(0.000799281 - model::transmission
		(params[0], 0., params[2], params[3], params[4], params[5])) < thresh);
	assert(abs(0.000801521 - model::static_conductance(params)) < thresh);
	assert(abs(0.000804088 - model::diff_conductance(params)) < thresh);

	params = {3., 1.4, -17., 0.67, 1.98, 1.3};
	assert(abs(0.00330201 - model::transmission
		(params[0], 0., params[2], params[3], params[4], params[5])) < thresh);
	assert(abs(0.00399841 - model::static_conductance(params)) < thresh);
	assert(abs(0.00480763 - model::diff_conductance(params)) < thresh);

	return 0;
}
