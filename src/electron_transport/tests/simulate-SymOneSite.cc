/* This file is a part of MolStat, which is distributed under the Creative
   Commons Attribution-NonCommercial 4.0 International Public License.
   MolStat (c) 2014, Northwestern University. */
/**
 * \file tests/simulate-SymOneSite.cc
 * \brief Test suite for the symmetric-coupling, single-site tight-binding
 *     model.
 *
 * \test Test suite for the symmetric-coupling, single-site tight-binding
 *     model.
 *
 * \author Matthew G.\ Reuter
 * \date October 2014
 */

#include <cstdio>
#include <cassert>
#include <cmath>
#include <array>
#include <map>

#include <general/random_distributions/constant.h>
#include <electron_transport/simulator_models/sym_one_site_simulate_model.h>

using namespace std;
using Model = molstat::transport::SymOneSiteSimulateModel;
constexpr size_t MPs{ 5 };

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
	const double thresh{ 1.e-6 };
	array<double, MPs> params;
	map<string, shared_ptr<molstat::RandomDistribution>> avail;

	// rig up some distributions for the constructor
	avail["ef"] = make_shared<molstat::ConstantDistribution>(0.);
	avail["v"] = make_shared<molstat::ConstantDistribution>(0.);
	avail["epsilon"] = make_shared<molstat::ConstantDistribution>(0.);
	avail["gamma"] = make_shared<molstat::ConstantDistribution>(0.);
	avail["a"] = make_shared<molstat::ConstantDistribution>(0.);

	// create a model
	Model model{ avail };

	// check known values for several parameter sets
	params[Model::Index_EF] = 0.;
	params[Model::Index_V] = 1.;
	params[Model::Index_epsilon] = -4.;
	params[Model::Index_gamma] = 0.8;
	params[Model::Index_a] = 0.;
	assert(abs(0.0384615 - model.transmission
		(params[Model::Index_EF], 0., params[Model::Index_epsilon],
		 params[Model::Index_gamma], params[Model::Index_a])) < thresh);
	assert(abs(0.0390172 - model.StaticG(params)) < thresh);
	assert(abs(0.0401438 - model.DiffG(params)) < thresh);
	assert(abs(params[Model::Index_V] - model.AppBias(params)) < thresh);

	params[Model::Index_EF] = 1.;
	params[Model::Index_V] = -0.4;
	params[Model::Index_epsilon] = -9.;
	params[Model::Index_gamma] = 0.4;
	params[Model::Index_a] = 0.;
	assert(abs(0.00159744 - model.transmission
		(params[Model::Index_EF], 0., params[Model::Index_epsilon],
		 params[Model::Index_gamma], params[Model::Index_a])) < thresh);
	assert(abs(0.00159808 - model.StaticG(params)) < thresh);
	assert(abs(0.00159936 - model.DiffG(params)) < thresh);
	assert(abs(params[Model::Index_V] - model.AppBias(params)) < thresh);

	params[Model::Index_EF] = 3.;
	params[Model::Index_V] = 1.4;
	params[Model::Index_epsilon] = -17.;
	params[Model::Index_gamma] = 0.67;
	params[Model::Index_a] = 0.;
	assert(abs(0.00112099 - model.transmission
		(params[Model::Index_EF], 0., params[Model::Index_epsilon],
		 params[Model::Index_gamma], params[Model::Index_a])) < thresh);
	assert(abs(0.00112236 - model.StaticG(params)) < thresh);
	assert(abs(0.00112511 - model.DiffG(params)) < thresh);
	assert(abs(params[Model::Index_V] - model.AppBias(params)) < thresh);

	params[Model::Index_EF] = 0.;
	params[Model::Index_V] = 1.;
	params[Model::Index_epsilon] = -4.;
	params[Model::Index_gamma] = 0.8;
	params[Model::Index_a] = -0.1;
	assert(abs(0.0384615 - model.transmission
		(params[Model::Index_EF], 0., params[Model::Index_epsilon],
		 params[Model::Index_gamma], params[Model::Index_a])) < thresh);
	assert(abs(0.0371825 - model.StaticG(params)) < thresh);
	assert(abs(0.0364382 - model.DiffG(params)) < thresh);
	assert(abs(params[Model::Index_V] - model.AppBias(params)) < thresh);

	params[Model::Index_EF] = 1.;
	params[Model::Index_V] = -0.4;
	params[Model::Index_epsilon] = -9.;
	params[Model::Index_gamma] = 0.4;
	params[Model::Index_a] = 1.;
	assert(abs(0.00159744 - model.transmission
		(params[Model::Index_EF], 0., params[Model::Index_epsilon],
		 params[Model::Index_gamma], params[Model::Index_a])) < thresh);
	assert(abs(0.00147765 - model.StaticG(params)) < thresh);
	assert(abs(0.00136520 - model.DiffG(params)) < thresh);
	assert(abs(params[Model::Index_V] - model.AppBias(params)) < thresh);

	params[Model::Index_EF] = 3.;
	params[Model::Index_V] = 1.4;
	params[Model::Index_epsilon] = -17.;
	params[Model::Index_gamma] = 0.67;
	params[Model::Index_a] = 0.24;
	assert(abs(0.00112099 - model.transmission
		(params[Model::Index_EF], 0., params[Model::Index_epsilon],
		 params[Model::Index_gamma], params[Model::Index_a])) < thresh);
	assert(abs(0.00116105 - model.StaticG(params)) < thresh);
	assert(abs(0.00120367 - model.DiffG(params)) < thresh);
	assert(abs(params[Model::Index_V] - model.AppBias(params)) < thresh);

	params[Model::Index_EF] = 3.;
	params[Model::Index_V] = 1.4;
	params[Model::Index_epsilon] = -17.;
	params[Model::Index_gamma] = 0.67;
	params[Model::Index_a] = -0.05;
	assert(abs(0.00112099 - model.transmission
		(params[Model::Index_EF], 0., params[Model::Index_epsilon],
		 params[Model::Index_gamma], params[Model::Index_a])) < thresh);
	assert(abs(0.00111445 - model.StaticG(params)) < thresh);
	assert(abs(0.00110948 - model.DiffG(params)) < thresh);
	assert(abs(params[Model::Index_V] - model.AppBias(params)) < thresh);

	return 0;
}
