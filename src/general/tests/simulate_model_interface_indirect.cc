/* This file is a part of MolStat, which is distributed under the Creative
   Commons Attribution-NonCommercial 4.0 International Public License.
   MolStat (c) 2014, Northwestern University. */
/**
 * \file simulate_model_interface_indirect.cc
 * \brief Test suite for the various MolStat classes that simulate data.
 *
 * \test Test suite for the various MolStat classes that simulate data. This
 *    test uses the indirect access functions.
 *
 * \author Matthew G.\ Reuter
 * \date October 2014
 */

#include <cassert>
#include "simulate_model_interface_observables.h"
#include "simulate_model_interface_models.h"
#include <general/random_distributions/rng.h>
#include <general/random_distributions/constant.h>
#include <general/simulator_tools/simulator.h>
#include <general/simulator_tools/simulate_model.h>

/**
 * \internal
 * \brief Main function for testing the various MolStat classes for
 *    simulating data. This test uses the indirect access functions.
 *
 * \param[in] argc The number of command-line arguments.
 * \param[in] argv The command-line arguments.
 * \return Exit status; 0 for normal.
 * \endinternal
 */
int main(int argc, char **argv) {
	constexpr double distvalue = 7.5;

	map<string, molstat::SimulateModelFactoryFunction> models;
	map<string, molstat::ObservableIndex> observables;

	molstat::gsl_rng_ptr r{ nullptr, &gsl_rng_free };

	// add the observables
	// must use emplace to avoid calling the allocator for type_index
	// (no default constructor)
	observables.emplace("obs1",
		molstat::GetObservableIndex<BasicObs1>() );
	observables.emplace("obs2",
		molstat::GetObservableIndex<BasicObs2>() );
	observables.emplace("obs3",
		molstat::GetObservableIndex<BasicObs3>() );
	observables.emplace("obs4",
		molstat::GetObservableIndex<BasicObs4>() );

	// add the model
	// use emplace for consistency
	models.emplace("basic",
		molstat::GetSimulateModelFactory<BasicTestModel>() );

	// get our factory
	molstat::SimulateModelFactory factory{ models.at("basic")() };

	// try to get the model. this should fail because we haven't specified
	// a distribution for the parameter "a"
	try{
		factory.getModel();

		assert(false);
	}
	catch(const molstat::MissingDistribution &e) {
		// should be here
	}

	// set the distribution
	// make sure the setDistribution function is case-insensitive
	factory.setDistribution("A",
		make_shared<molstat::ConstantDistribution>(distvalue));

	// wrap the model into our simulator
	molstat::Simulator sim( factory.getModel() );

	// try to simulate data. this should fail because we haven't set any
	// observables or distributions yet.
	try {
		sim.simulate(r);

		assert(false);
	}
	catch(const molstat::NoObservables &e) {
		// we should be here
	}

	// try to set an observable with a bad index.
	try {
		// only 0 is available right now...
		sim.setObservable(1, observables.at("obs1"));

		assert(false);
	}
	catch(const out_of_range &e) {
		// should be here
	}

	// try to set Observable4
	// this should fail because TestModel doesn't implement Observable4
	try {
		sim.setObservable(0, observables.at("obs4"));

		assert(false);
	}
	catch(const molstat::IncompatibleObservable &e) {
		// should be here
	}

	// now, let's actually set an observable
	sim.setObservable(0, observables.at("obs1"));

	valarray<double> data = sim.simulate(r);
	assert(abs(data[0] - distvalue) < 1.e-6);

	// set another observable
	sim.setObservable(1, observables.at("obs2"));

	// verify the set of observables generated...
	data = sim.simulate(r);
	assert(abs(data[0] - distvalue) < 1.e-6);
	assert(abs(data[1] - BasicTestModel::obs2value) < 1.e-6);

	// change an observable and recheck
	sim.setObservable(0, observables.at("obs2"));
	data = sim.simulate(r);
	assert(abs(data[0] - BasicTestModel::obs2value) < 1.e-6);
	assert(abs(data[1] - BasicTestModel::obs2value) < 1.e-6);

	// now load in Observable3
	sim.setObservable(1, observables.at("obs3"));

	try {
		// Observable3 should throw an exception; make sure we catch it
		sim.simulate(r);
		assert(false);
	}
	catch(int i) {
		assert(i == BasicTestModel::obs3except);
	}

	return 0;
}
