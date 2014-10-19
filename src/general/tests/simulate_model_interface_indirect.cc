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
#include "simulate_model_interface_classes.h"
#include <general/random_distributions/rng.h>
#include <general/random_distributions/constant.h>
#include <general/simulator_tools/simulator.h>

/**
 * \internal
 * \brief Main function for testing the molstat::SimulateObservables,
 *    molstat::SimulateModel, and molstat::Observable templates.
 *
 * \param[in] argc The number of command-line arguments.
 * \param[in] argv The command-line arguments.
 * \return Exit status; 0 for normal.
 * \endinternal
 */
int main(int argc, char **argv) {
	map<string, molstat::SimulateModelFactory> models;
	map<string, molstat::ObservableIndex> observables;

	molstat::gsl_rng_ptr r{ nullptr, &gsl_rng_free };

	// add the observables
	// must use emplace to avoid calling the allocator for type_index
	// (no default constructor)
	observables.emplace(make_pair("obs1",
		molstat::GetObservableIndex<Observable1>()) );
	observables.emplace(make_pair("obs2",
		molstat::GetObservableIndex<Observable2>()) );
	observables.emplace(make_pair("obs3",
		molstat::GetObservableIndex<Observable3>()) );
	observables.emplace(make_pair("obs4",
		molstat::GetObservableIndex<Observable4>()) );

	// add the model
	// use emplace for consistency
	models.emplace(make_pair("test",
		molstat::GetSimulateModelFactory<TestModel>() ));

	// make our model using the functor
	shared_ptr<molstat::SimulateModel> model { models.at("test")() };

	// wrap the model into our simulator
	molstat::Simulator sim { model };

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

	// try to simulate data. this should now fail because we haven't specified
	// a distribution for the parameter "a"
	try {
		sim.simulate(r);

		assert(false);
	}
	catch(const molstat::MissingDistribution &e) {
		// should be here
	}

	// set the distribution
	// make sure the setDistribution function is case-insensitive
	model->setDistribution("A",
		make_shared<molstat::ConstantDistribution>(distvalue));

	// verify the set
	valarray<double> data = sim.simulate(r);
	assert(abs(data[0] - distvalue) < 1.e-6);

	model->setDistribution("a",
		make_shared<molstat::ConstantDistribution>(distvalue));

	data = sim.simulate(r);
	assert(abs(data[0] - distvalue) < 1.e-6);

	// set another observable
	sim.setObservable(1, observables.at("obs2"));

	// verify the set of observables generated...
	data = sim.simulate(r);
	assert(abs(data[0] - distvalue) < 1.e-6);
	assert(abs(data[1] - constvalue) < 1.e-6);

	// change an observable and recheck
	sim.setObservable(0, observables.at("obs2"));
	data = sim.simulate(r);
	assert(abs(data[0] - constvalue) < 1.e-6);
	assert(abs(data[1] - constvalue) < 1.e-6);

	// now load in Observable3
	sim.setObservable(1, observables.at("obs3"));

	try {
		// Observable3 should throw an exception; make sure we catch it
		sim.simulate(r);
		assert(false);
	}
	catch(int i) {
		assert(i == exceptvalue);
	}

	return 0;
}
