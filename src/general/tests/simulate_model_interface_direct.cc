/* This file is a part of MolStat, which is distributed under the Creative
   Commons Attribution-NonCommercial 4.0 International Public License.
   MolStat (c) 2014, Northwestern University. */
/**
 * \file simulate_model_interface_direct.cc
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
	unique_ptr<molstat::Simulator<3>> sim;
	map<string, shared_ptr<molstat::RandomDistribution>> parameters;
	molstat::gsl_rng_ptr r(nullptr, &gsl_rng_free);

	try {
		sim = molstat::Simulator<3>::Factory<TestModel>(parameters);

		// shouldn't be here because we didn't specify the parameter "a"
		assert(false);
	}
	catch(const runtime_error &e) {
		// we should be here!
	}

	parameters["a"] = make_shared<molstat::ConstantDistribution>(distvalue);

	try {
		// this time it should work!
		sim = molstat::Simulator<3>::Factory<TestModel>(parameters);
	}
	catch(const runtime_error &e) {
		assert(false);
	}

	// try to set an observable for Observable4.
	// This should fail (TestModel doesn't implement Observable4)
	try {
		sim->setObservable<Observable4>(0);
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
		sim->setObservable<Observable1>(3);
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
		sim->setObservable<Observable1>(0);
		sim->setObservable<Observable2>(2);
	}
	catch(const runtime_error &e) {
		// this should have worked...
		assert(false);
	}
	catch(const out_of_range &e) {
		assert(false);
	}

	// verify the set of observables generated...
	array<double, 3> data = sim->simulate(r);
	assert(abs(data[0] - distvalue) < 1.e-6);
	assert(abs(data[1] - 0.) < 1.e-6);
	assert(abs(data[2] - constvalue) < 1.e-6);

	// now load in Observable3
	try {
		sim->setObservable<Observable3>(1);
	}
	catch(const runtime_error &e) {
		// the set should have worked.
		assert(false);
	}
	catch(const out_of_range &e) {
		assert(false);
	}

	try {
		// Observable 3 should throw and exception; make sure we catch it
		data = sim->simulate(r);
		assert(false);
	}
	catch(int i) {
		assert(i == exceptvalue);
	}

	// check the ordering of distributions within a model
	parameters["b"] = make_shared<molstat::ConstantDistribution>(distvalue+1.);
	parameters["c"] = make_shared<molstat::ConstantDistribution>(distvalue-1.);

	try {
		// attempt to create the simulator; this should fail because item 2 is
		// missing in the construction of a FailedMapModel
		sim = molstat::Simulator<3>::Factory<FailedMapModel>(parameters);

		assert(false);
	}
	catch(const out_of_range &e) {
		// should be here
	}
	catch(const runtime_error &e) {
		assert(false);
	}

	try {
		unique_ptr<GoodMapModel> model{ new GoodMapModel(parameters) };

		// make sure the correct distributions are in the correct places
		assert(model->order[0] == "a");
		assert(model->order[1] == "c");
		assert(model->order[2] == "a");
	}
	catch(...) {
		assert(false);
	}

	return 0;
}
