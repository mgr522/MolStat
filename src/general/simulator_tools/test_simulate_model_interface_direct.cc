/* This file is a part of MolStat, which is distributed under the Creative
   Commons Attribution-NonCommercial 4.0 International Public License.
   MolStat (c) 2014, Northwestern University. */
/**
 * \file test_simulate_model_interface_direct.cc
 * \brief Test suite for the molstat::SimulateObservables,
 *    molstat::SimulateModel, and molstat::Observable templates.
 *
 * \test Tests the molstat::SimulateObservables, molstat::SimulateModel,
 *    and molstat::Observable templates using direct access to the
 *    molstat::SimulatorFactory object.
 *
 * \author Matthew G.\ Reuter
 * \date October 2014
 */

#include "test_simulate_model_interface_classes.h"
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
	molstat::SimulatorFactory<3> factory;
	shared_ptr<molstat::Simulator<3>> sim;
	map<string, shared_ptr<molstat::RandomDistribution>> parameters;
	molstat::gsl_rng_ptr r(nullptr, &gsl_rng_free);

	try {
		factory = molstat::SimulatorFactory<3>
			::makeFactory<TestModel>(parameters);

		// shouldn't be here because we didn't specify the parameter "a"
		assert(false);
	}
	catch(const runtime_error &e) {
		// we should be here!
	}

	parameters["a"] = make_shared<molstat::ConstantDistribution>(distvalue);

	try {
		// this time it should work!
		factory = molstat::SimulatorFactory<3>
			::makeFactory<TestModel>(parameters);
	}
	catch(const runtime_error &e) {
		assert(false);
	}

	// try to set an observable for Observable4.
	// This should fail (TestModel doesn't implement Observable4)
	try {
		factory.setObservable<Observable4>(0);
		assert(false);
	}
	catch(const runtime_error &e) {
		// should be here
	}
	catch(const out_of_range &e) {
		assert(false);
	}

	// try to set an observable to a bad index
	try {
		factory.setObservable<Observable1>(3);
		assert(false);
	}
	catch(const out_of_range &e) {
		// should be here
	}
	catch(const runtime_error &e) {
		assert(false);
	}

	// now set observables for Observable1 and Observable2.
	try {
		factory.setObservable<Observable1>(0);
		factory.setObservable<Observable2>(2);
	}
	catch(const runtime_error &e) {
		// this should have worked...
		assert(false);
	}
	catch(const out_of_range &e) {
		assert(false);
	}

	// cast to the simulator now that setup is complete
	sim = factory.create();
	assert(sim != nullptr);

	// verify the set of observables generated...
	array<double, 3> data = sim->simulate(r);
	assert(abs(data[0] - distvalue) < 1.e-6);
	assert(abs(data[1] - 0.) < 1.e-6);
	assert(abs(data[2] - constvalue) < 1.e-6);

	// create a new simulator that uses Observable3
	try {
		// this time it should work!
		factory = molstat::SimulatorFactory<3>
			::makeFactory<TestModel>(parameters);
	}
	catch(const runtime_error &e) {
		assert(false);
	}

	// now load in Observable3
	try {
		factory.setObservable<Observable1>(0);
		factory.setObservable<Observable3>(1);
		factory.setObservable<Observable2>(2);
	}
	catch(const runtime_error &e) {
		// the set should have worked.
		assert(false);
	}
	catch(const out_of_range &e) {
		assert(false);
	}

	sim = factory.create();
	assert(sim != nullptr);

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
