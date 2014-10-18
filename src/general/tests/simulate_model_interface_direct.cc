/* This file is a part of MolStat, which is distributed under the Creative
   Commons Attribution-NonCommercial 4.0 International Public License.
   MolStat (c) 2014, Northwestern University. */
/**
 * \file simulate_model_interface_direct.cc
 * \brief Test suite for the various MolStat classes for simulating data.
 *
 * \test Test suite for the various MolStat classes for simulating data.
 *
 * \author Matthew G.\ Reuter
 * \date October 2014
 */

#include <cassert>
#include <typeinfo>
#include <typeindex>
#include "simulate_model_interface_classes.h"
#include <general/random_distributions/rng.h>
#include <general/random_distributions/constant.h>
#include <general/simulator_tools/simulator.h>

/**
 * \internal
 * \brief Main function for testing the various MolStat classes for
 *    simulating data.
 *
 * \param[in] argc The number of command-line arguments.
 * \param[in] argv The command-line arguments.
 * \return Exit status; 0 for normal.
 * \endinternal
 */
int main(int argc, char **argv) {
	molstat::gsl_rng_ptr r{ nullptr, &gsl_rng_free }; // dummy rng
	shared_ptr<molstat::SimulateModel> model { make_shared<TestModel>() };

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
		sim.setObservable(1, type_index{ typeid(Observable1) });

		assert(false);
	}
	catch(const out_of_range &e) {
		// should be here
	}

	// try to set Observable4
	// this should fail because TestModel doesn't implement Observable4
	try {
		sim.setObservable(0, type_index{ typeid(Observable4) });

		assert(false);
	}
	catch(const molstat::IncompatibleObservable &e) {
		// should be here
	}

	// now, let's actually set an observable
	sim.setObservable(0, type_index{ typeid(Observable1) });

	// try to simulate data. this should now fail because we haven't specified
	// a distribution for the parameter "a"
	try{
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
	sim.setObservable(1, type_index{ typeid(Observable2) });

	// verify the set of observables generated...
	data = sim.simulate(r);
	assert(abs(data[0] - distvalue) < 1.e-6);
	assert(abs(data[1] - constvalue) < 1.e-6);

	// change an observable and recheck
	sim.setObservable(0, type_index{ typeid(Observable2) });
	data = sim.simulate(r);
	assert(abs(data[0] - constvalue) < 1.e-6);
	assert(abs(data[1] - constvalue) < 1.e-6);

	// now load in Observable3
	sim.setObservable(1, type_index{ typeid(Observable3) });

	try {
		// Observable 3 should throw and exception; make sure we catch it
		sim.simulate(r);
		assert(false);
	}
	catch(int i) {
		assert(i == exceptvalue);
	}

	return 0;
}
