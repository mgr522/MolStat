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
#include "simulate_model_interface_observables.h"
#include "simulate_model_interface_models.h"
#include <general/random_distributions/rng.h>
#include <general/random_distributions/constant.h>
#include <general/simulator_tools/simulator.h>
#include <general/simulator_tools/simulate_model.h>

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
	molstat::SimulateModelFactory factory 
		{ molstat::SimulateModelFactory::makeFactory<BasicTestModel>() };

	molstat::gsl_rng_ptr r{ nullptr, &gsl_rng_free }; // dummy rng

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
		sim.setObservable(1, type_index{ typeid(BasicObs1) });

		assert(false);
	}
	catch(const out_of_range &e) {
		// should be here
	}

	// try to set BasicObs4
	// this should fail because TestModel doesn't implement BasicObs4
	try {
		sim.setObservable(0, type_index{ typeid(BasicObs4) });

		assert(false);
	}
	catch(const molstat::IncompatibleObservable &e) {
		// should be here
	}

	// now, let's actually set an observable
	sim.setObservable(0, type_index{ typeid(BasicObs1) });

	valarray<double> data = sim.simulate(r);
	assert(abs(data[0] - distvalue) < 1.e-6);

	// set another observable
	sim.setObservable(1, type_index{ typeid(BasicObs2) });

	// verify the set of observables generated...
	data = sim.simulate(r);
	assert(abs(data[0] - distvalue) < 1.e-6);
	assert(abs(data[1] - constvalue) < 1.e-6);

	// change an observable and recheck
	sim.setObservable(0, type_index{ typeid(BasicObs2) });
	data = sim.simulate(r);
	assert(abs(data[0] - constvalue) < 1.e-6);
	assert(abs(data[1] - constvalue) < 1.e-6);

	// now load in Observable3
	sim.setObservable(1, type_index{ typeid(BasicObs3) });

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
