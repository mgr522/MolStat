/* This file is a part of MolStat, which is distributed under the Creative
   Commons Attribution-NonCommercial 4.0 International Public License.
   MolStat (c) 2014, Northwestern University. */
/**
 * \file simulate_model_interface_indirect.cc
 * \brief Test suite for the molstat::SimulateObservables,
 *    molstat::SimulateModel, and molstat::Observable templates.
 *
 * \test Tests the molstat::SimulateObservables, molstat::SimulateModel,
 *    and molstat::Observable templates using the access functions.
 *
 * \author Matthew G.\ Reuter
 * \date October 2014
 */

#include "simulate_model_interface_classes.h"
#include <general/random_distributions/constant.h>

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
	map<string, molstat::SimulatorFactory<3>> models;
	map<string, molstat::ObservableSetter<3>> observables;

	shared_ptr<molstat::Simulator<3>> sim;
	map<string, shared_ptr<molstat::RandomDistribution>> parameters;
	molstat::gsl_rng_ptr r(nullptr, &gsl_rng_free);
	bool flag;

	// add the model
	models["test"] = molstat::GetSimulatorFactory<3, TestModel>();

	// add the observables
	observables["obs1"] = molstat::GetObservableSetter<3, Observable1>();
	observables["obs2"] = molstat::GetObservableSetter<3, Observable2>();
	observables["obs3"] = molstat::GetObservableSetter<3, Observable3>();
	observables["obs4"] = molstat::GetObservableSetter<3, Observable4>();

	try {
		// use a flag to make sure out_of_range is thrown at the correct spot
		flag = false;
		auto factory = models.at("test");
		flag = true;
		sim = factory(parameters);

		// shouldn't be here because we didn't specify the parameter "a"
		assert(false);
	}
	catch(const out_of_range &e) {
		// we should be here, but from the factory(parameters) call
		assert(flag);
	}

	// create an entry in parameters["a"], but leave it set to nullptr
	parameters["a"];
	try {
		auto factory = models.at("test");
		sim = factory(parameters);

		// should throw because "a" associates with nullptr
		assert(false);
	}
	catch(const runtime_error &e) {
		// should be here
	}

	parameters["a"] = make_shared<molstat::ConstantDistribution>(distvalue);

	// this time it should work!
	sim = models.at("test")(parameters);

	// try to set an observable for Observable4.
	// This should fail (TestModel doesn't implement Observable4)
	try {
		observables.at("obs4")(sim.get(), 0);
		assert(false);
	}
	catch(const runtime_error &e) {
		// should be here
	}
	
	// try to set an observable to a bad index
	// need a flag to make sure we get out_of_range for the correct reason
	try {
		flag = false;
		auto functor = observables.at("obs1");
		flag = true;
		functor(sim.get(), 3);
		assert(false);
	}
	catch(const out_of_range &e) {
		// should be here
		assert(flag);
	}

	// now set observables for Observable1 and Observable2.
	observables.at("obs1")(sim.get(), 0);
	observables.at("obs2")(sim.get(), 2);

	// verify the set of observables generated...
	array<double, 3> data = sim->simulate(r);
	assert(abs(data[0] - distvalue) < 1.e-6);
	assert(abs(data[1] - 0.) < 1.e-6);
	assert(abs(data[2] - constvalue) < 1.e-6);

	// now load in Observable3
	observables.at("obs3")(sim.get(), 1);

	try {
		// Observable 3 should throw and exception; make sure we catch it
		data = sim->simulate(r);
		assert(false);
	}
	catch(int i) {
		assert(i == exceptvalue);
	}

	return 0;
}