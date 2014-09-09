/* This file is a part of MolStat, which is distributed under the Creative
   Commons Attribution-NonCommercial 4.0 International Public License.
   MolStat (c) 2014, Northwestern University. */
/**
 * \file tests/symmetric_one_site.cc
 * \brief Test suite for the symmetric-coupling, single-site tight-binding
 *     model.
 *
 * \test Test suite for the symmetric-coupling, single-site tight-binding
 *     model.
 *
 * \author Matthew G.\ Reuter
 * \date September 2014
 */

#include <cstdio>
#include <assert.h>
#include <cmath>
#include <vector>

#include "../simulator_models/sym_one_site_simulate_model.h"

using namespace std;

/**
 * \internal
 * \brief Main function for testing the symmetric-coupling, one-site model.
 *
 * \param[in] argc The number of command-line arguments.
 * \param[in] argv The command-line arguments.
 * \return Exit status: 0 if the code passes the test, non-zero otherwise.
 * \endinternal
 */
int main(int argc, char **argv) {
	const double thresh = 1.0e-6;
	typedef SymOneSiteSimulateModel model;
	vector<double> params(5);

	// params[0] = EF, [1] = V, [2] = eps, [3] = gamma, [4] = a

	// check known values for several parameter sets
	params = {0., 1., -4., 0.8, 0.};
	assert(abs(0.0384615 - model::transmission
		(params[0], 0., params[2], params[3], params[4])) < thresh);
	assert(abs(0.0390172 - model::static_conductance(params)) < thresh);
	assert(abs(0.0401438 - model::diff_conductance(params)) < thresh);

	params = {1., -0.4, -9., 0.4, 0.};
	assert(abs(0.00159744 - model::transmission
		(params[0], 0., params[2], params[3], params[4])) < thresh);
	assert(abs(0.00159808 - model::static_conductance(params)) < thresh);
	assert(abs(0.00159936 - model::diff_conductance(params)) < thresh);

	params = {3., 1.4, -17., 0.67, 0.};
	assert(abs(0.00112099 - model::transmission
		(params[0], 0., params[2], params[3], params[4])) < thresh);
	assert(abs(0.00112236 - model::static_conductance(params)) < thresh);
	assert(abs(0.00112511 - model::diff_conductance(params)) < thresh);

	params = {0., 1., -4., 0.8, -0.1};
	assert(abs(0.0384615 - model::transmission
		(params[0], 0., params[2], params[3], params[4])) < thresh);
	assert(abs(0.0371825 - model::static_conductance(params)) < thresh);
	assert(abs(0.0364382 - model::diff_conductance(params)) < thresh);

	params = {1., -0.4, -9., 0.4, 1.};
	assert(abs(0.00159744 - model::transmission
		(params[0], 0., params[2], params[3], params[4])) < thresh);
	assert(abs(0.00147765 - model::static_conductance(params)) < thresh);
	assert(abs(0.00136520 - model::diff_conductance(params)) < thresh);

	params = {3., 1.4, -17., 0.67, 0.24};
	assert(abs(0.00112099 - model::transmission
		(params[0], 0., params[2], params[3], params[4])) < thresh);
	assert(abs(0.00116105 - model::static_conductance(params)) < thresh);
	assert(abs(0.00120367 - model::diff_conductance(params)) < thresh);

	params = {3., 1.4, -17., 0.67, -0.05};
	assert(abs(0.00112099 - model::transmission
		(params[0], 0., params[2], params[3], params[4])) < thresh);
	assert(abs(0.00111455 - model::static_conductance(params)) < thresh);
	assert(abs(0.00110948 - model::diff_conductance(params)) < thresh);

	return 0;
}
