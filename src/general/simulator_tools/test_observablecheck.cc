/* This file is a part of MolStat, which is distributed under the Creative
   Commons Attribution-NonCommercial 4.0 International Public License.
   MolStat (c) 2014, Northwestern University. */
/**
 * \file test_observablecheck.cc
 * \brief Test suite for the ::ObservableCheck template.
 *
 * \test Tests the ::ObservableCheck template.
 *
 * \author Matthew G.\ Reuter
 * \date September 2014
 */

#include <memory>
#include <cstdio>
#include <map>
#include <string>
#include <functional>
#include <array>
#include <assert.h>

#include <electron_transport/simulator_models/transport_observables.h>
#include <general/random_distributions/constant.h>
#include "simulate_model_interface.h"

using namespace std;

/**
 * \internal
 * \brief Dummy class for testing ::ObservableCheck.
 * \endinternal
 */
class ObservableCheckClass : public molstat::SimulateModel,
	public molstat::StaticConductance {
public:
	ObservableCheckClass() = delete;
	virtual ~ObservableCheckClass() = default;

	ObservableCheckClass(
		const map<string, shared_ptr<molstat::RandomDistribution>> &avail)
		: molstat::SimulateModel(avail, { "a" }) {}

	virtual array<double, 2> StaticG(molstat::gsl_rng_ptr &r) const override {
		throw 4;
		return { {0., 0.} };
	}
};

/**
 * \internal
 * \brief Main function for testing the ::ObservableCheck template.
 *
 * \param[in] argc The number of command-line arguments.
 * \param[in] argv The command-line arguments.
 * \return Exit status; 0 for normal.
 * \endinternal
 */
int main(int argc, char **argv) {
	shared_ptr<molstat::SimulateModel> model;
	molstat::Observable<2> obs;
	function<molstat::Observable<2>(shared_ptr<molstat::SimulateModel>) > func;
	map<string, shared_ptr<molstat::RandomDistribution>> parameters;

	try {
		model = make_shared<ObservableCheckClass>(parameters);

		// shouldn't be here because we didn't specify the parameter "a"
		assert(false);
	}
	catch(const runtime_error &e) {
		// we should be here!
	}

	parameters["a"] = make_shared<molstat::ConstantDistribution>(0.);

	try {
		// this time it should work!
		model = make_shared<ObservableCheckClass>(parameters);
	}
	catch(const runtime_error &e) {
		assert(false);
	}

	// try to check for DifferentialConductance. This should fail.
	try {
		func = molstat::ObservableCheck(
			&molstat::DifferentialConductance::DiffG);
		func(model); // this should throw because the model doesn't have DiffG
		assert(false);
	}
	catch(const runtime_error &e) {
		// should be here
	}

	// now check for StaticConductance.
	try {
		func = molstat::ObservableCheck(
			&molstat::StaticConductance::StaticG);
		obs = func(model);
	}
	catch(const runtime_error &e) {
		// this should have worked...
		assert(false);
	}

	// verify that obs is pointing to the correct function
	// note that this test is not conclusive...
	try {
		molstat::gsl_rng_ptr r(nullptr, &gsl_rng_free);
		obs(r); // this should throw 4
		assert(false);
	}
	catch(int i) {
		assert(i == 4);
	}

	return 0;
}
