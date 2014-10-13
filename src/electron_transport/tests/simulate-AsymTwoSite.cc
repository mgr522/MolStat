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
 * \date October 2014
 */

#include <cstdio>
#include <cassert>
#include <cmath>
#include <array>
#include <map>

#include <general/random_distributions/constant.h>
#include <electron_transport/simulator_models/asym_two_site_simulate_model.h>

using namespace std;
using Model = molstat::transport::AsymTwoSiteSimulateModel;
constexpr size_t MPs{ 6 };

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
	array<double, MPs> params;
	map<string, shared_ptr<molstat::RandomDistribution>> avail;

	// rig up some distributions for the constructor
	avail["ef"] = make_shared<molstat::ConstantDistribution>(0.);
	avail["v"] = make_shared<molstat::ConstantDistribution>(0.);
	avail["epsilon"] = make_shared<molstat::ConstantDistribution>(0.);
	avail["gammal"] = make_shared<molstat::ConstantDistribution>(0.);
	avail["gammar"] = make_shared<molstat::ConstantDistribution>(0.);
	avail["beta"] = make_shared<molstat::ConstantDistribution>(0.);

	// create a model
	Model model{ avail };

	// check known values for several parameter sets
	params[Model::Index_EF] = 0.;
	params[Model::Index_V] = 1.;
	params[Model::Index_epsilon] = -4.;
	params[Model::Index_gammaL] = 0.8;
	params[Model::Index_gammaR] = 1.;
	params[Model::Index_beta] = -3.;
	assert(abs(0.121622 - model.transmission
		(params[Model::Index_EF], 0., params[Model::Index_epsilon],
		 params[Model::Index_gammaL], params[Model::Index_gammaR],
		 params[Model::Index_beta])) < thresh);
	assert(abs(0.149936 - model.StaticG(params)) < thresh);
	assert(abs(0.213248 - model.DiffG(params)) < thresh);
	assert(abs(params[Model::Index_V] - model.AppBias(params)) < thresh);

	params[Model::Index_EF] = 1.;
	params[Model::Index_V] = -0.4;
	params[Model::Index_epsilon] = -3.;
	params[Model::Index_gammaL] = 0.4;
	params[Model::Index_gammaR] = 0.2;
	params[Model::Index_beta] = -0.8;
	assert(abs(0.000216257 - model.transmission
		(params[Model::Index_EF], 0., params[Model::Index_epsilon],
		 params[Model::Index_gammaL], params[Model::Index_gammaR],
		 params[Model::Index_beta])) < thresh);
	assert(abs(0.000218231 - model.StaticG(params)) < thresh);
	assert(abs(0.000222203 - model.DiffG(params)) < thresh);
	assert(abs(params[Model::Index_V] - model.AppBias(params)) < thresh);

	params[Model::Index_EF] = -1.;
	params[Model::Index_V] = 1.4;
	params[Model::Index_epsilon] = 5.;
	params[Model::Index_gammaL] = 0.67;
	params[Model::Index_gammaR] = 1.98;
	params[Model::Index_beta] = -1.6;
	assert(abs(0.00292927 - model.transmission
		(params[Model::Index_EF], 0., params[Model::Index_epsilon],
		 params[Model::Index_gammaL], params[Model::Index_gammaR],
		 params[Model::Index_beta])) < thresh);
	assert(abs(0.00308371 - model.StaticG(params)) < thresh);
	assert(abs(0.00340305 - model.DiffG(params)) < thresh);
	assert(abs(params[Model::Index_V] - model.AppBias(params)) < thresh);

	return 0;
}
