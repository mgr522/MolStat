/* This file is a part of MolStat, which is distributed under the Creative
   Commons Attribution-NonCommercial 4.0 International Public License.

   (c) 2014 Northwestern University. */

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
 * \brief Main function for testing the various MolStat classes for
 *    simulating data.
 *
 * \param[in] argc The number of command-line arguments.
 * \param[in] argv The command-line arguments.
 * \return Exit status; 0 for normal.
 */
int main(int argc, char **argv)
{
	constexpr double distvalue1 = 7.5;
	constexpr double distvalue2 = -6.3;
	constexpr double distvalue3 = 4.1;
	constexpr double distvalue4 = 1.2;
	constexpr double distvalue5 = -14.1;
	constexpr double distvalue6 = 2.9;

	molstat::SimulateModelFactory basic_factory 
		{ molstat::SimulateModelFactory::makeFactory<BasicTestModel>() };

	molstat::Engine engine; // dummy random number engine

	// try to get the model. this should fail because we haven't specified
	// a distribution for the parameter "a"
	try
	{
		basic_factory.getModel();
		assert(false);
	}
	catch(const molstat::MissingDistribution &e)
	{
		// should be here
	}

	// set the distribution
	// make sure the setDistribution function is case-insensitive
	basic_factory.setDistribution("A",
		make_shared<molstat::ConstantDistribution>(distvalue1));

	// wrap the model into our simulator
	// save the model independently for use in a later test
	shared_ptr<molstat::SimulateModel> basic_model{ basic_factory.getModel() };
	molstat::Simulator basic_sim{ basic_model };

	// try to simulate data. this should fail because we haven't set any
	// observables or distributions yet.
	try
	{
		basic_sim.simulate(engine);
		assert(false);
	}
	catch(const molstat::NoObservables &e)
	{
		// we should be here
	}

	// try to set an observable with a bad index.
	try
	{
		// only 0 is available right now...
		basic_sim.setObservable(1, type_index{ typeid(BasicObs1) });
		assert(false);
	}
	catch(const out_of_range &e)
	{
		// should be here
	}

	// try to set BasicObs4
	// this should fail because TestModel doesn't implement BasicObs4
	try
	{
		basic_sim.setObservable(0, type_index{ typeid(BasicObs4) });
		assert(false);
	}
	catch(const molstat::IncompatibleObservable &e)
	{
		// should be here
	}

	// now, let's actually set an observable
	basic_sim.setObservable(0, type_index{ typeid(BasicObs1) });

	valarray<double> data = basic_sim.simulate(engine);
	assert(abs(data[0] - distvalue1) < 1.e-6);

	// set another observable
	basic_sim.setObservable(1, type_index{ typeid(BasicObs2) });

	// verify the set of observables generated...
	data = basic_sim.simulate(engine);
	assert(abs(data[0] - distvalue1) < 1.e-6);
	assert(abs(data[1] - BasicTestModel::obs2value) < 1.e-6);

	// change an observable and recheck
	basic_sim.setObservable(0, type_index{ typeid(BasicObs2) });
	data = basic_sim.simulate(engine);
	assert(abs(data[0] - BasicTestModel::obs2value) < 1.e-6);
	assert(abs(data[1] - BasicTestModel::obs2value) < 1.e-6);

	// now load in Observable3
	basic_sim.setObservable(1, type_index{ typeid(BasicObs3) });

	try
	{
		// Observable3 should throw an exception; make sure we catch it
		basic_sim.simulate(engine);
		assert(false);
	}
	catch(int i)
	{
		assert(i == BasicTestModel::obs3except);
	}

	// now start testing functionality of composite models.

	// as is, the factory object is still set to the basic model.
	// try to set a submodel.
	try
	{
		basic_factory.addSubmodel(nullptr);
		assert(false);
	}
	catch(const molstat::NotCompositeSimulateModel &e)
	{
		// should be here
	}

	// make factories for the composite models
	molstat::SimulateModelFactory cfactory_add
		{ molstat::SimulateModelFactory::makeFactory<CompositeTestModelAdd>() };
	molstat::SimulateModelFactory cfactory_mult
		{ molstat::SimulateModelFactory::
			makeFactory<CompositeTestModelMultiply>() };

	// try to set an invalid submodel (wrong type)
	try
	{
		cfactory_add.addSubmodel(basic_model);
		assert(false);
	}
	catch(const molstat::IncompatibleSubmodel &e)
	{
		// should be here
	}

	// specify the required distributions
	cfactory_add.setDistribution("ef",
		make_shared<molstat::ConstantDistribution>(distvalue1));
	cfactory_add.setDistribution("v",
		make_shared<molstat::ConstantDistribution>(distvalue2));
	cfactory_mult.setDistribution("ef",
		make_shared<molstat::ConstantDistribution>(distvalue1));
	cfactory_mult.setDistribution("v",
		make_shared<molstat::ConstantDistribution>(distvalue2));

	// try to get access to the model
	// this should fail because we haven't specified any submodels yet
	try
	{
		cfactory_add.getModel();
		assert(false);
	}
	catch(const molstat::NoSubmodels &e)
	{
		// should be here
	}

	// create submodels
	molstat::SimulateModelFactory subfactory1
		{ molstat::SimulateModelFactory::makeFactory<CompositeSubModel>() };
	molstat::SimulateModelFactory subfactory2
		{ molstat::SimulateModelFactory::makeFactory<CompositeSubModel>() };

	// add a distribution to the submodels
	subfactory1.setDistribution("eps",
		make_shared<molstat::ConstantDistribution>(distvalue3));
	subfactory2.setDistribution("gamma",
		make_shared<molstat::ConstantDistribution>(distvalue6));

	// try to get a model
	// this should still fail because one of the distributions is unspecified
	try
	{
		subfactory1.getModel();
		assert(false);
	}
	catch(const molstat::MissingDistribution &e)
	{
		// should be here
	}

	try
	{
		subfactory2.getModel();
		assert(false);
	}
	catch(const molstat::MissingDistribution &e)
	{
		// should be here
	}

	subfactory1.setDistribution("gamma",
		make_shared<molstat::ConstantDistribution>(distvalue4));
	subfactory2.setDistribution("eps",
		make_shared<molstat::ConstantDistribution>(distvalue5));

	// get the submodels
	shared_ptr<molstat::SimulateModel> submodel1 { subfactory1.getModel() };
	shared_ptr<molstat::SimulateModel> submodel2 { subfactory2.getModel() };

	// add one submodel to the additive model
	cfactory_add.addSubmodel(submodel1);

	// test the composite observable with just one submodel
	{
		// make a simulator
		molstat::Simulator sim_add{ cfactory_add.getModel() };
		// add the observables
		sim_add.setObservable(0, type_index{ typeid(BasicObs4) });
		sim_add.setObservable(1, type_index{ typeid(BasicObs1) });
		data = sim_add.simulate(engine);
		assert(abs(data[0] - (distvalue1 + distvalue2)) < 1.e-6);
		assert(abs(data[1] - (distvalue2 * (distvalue3 - distvalue4))) < 1.e-6);
	}

	// add the rest of the submodels
	cfactory_add.addSubmodel(submodel2);
	cfactory_mult.addSubmodel(submodel1);
	cfactory_mult.addSubmodel(submodel2);

	// make simulators for the two composite models
	molstat::Simulator sim_add{ cfactory_add.getModel() };
	molstat::Simulator sim_mult{ cfactory_mult.getModel() };

	// set observables
	// first try one that fails
	try
	{
		sim_add.setObservable(0, type_index{ typeid(BasicObs2) });
		assert(false);
	}
	catch(const molstat::IncompatibleObservable &e)
	{
		// should be here
	}

	// now for observables that work
	sim_add.setObservable(0, type_index{ typeid(BasicObs4) });
	sim_add.setObservable(1, type_index{ typeid(BasicObs1) });
	sim_mult.setObservable(0, type_index{ typeid(BasicObs4) });
	sim_mult.setObservable(1, type_index{ typeid(BasicObs1) });

	// simulate data to test the combination of values
	data = sim_add.simulate(engine);
	assert(abs(data[0] - (distvalue1 + distvalue2)) < 1.e-6);
	assert(abs(data[1] - (distvalue2 * (distvalue3 - distvalue4) + // sum
	                      distvalue2 * (distvalue5 - distvalue6))) < 1.e-6);
	data = sim_mult.simulate(engine);
	assert(abs(data[0] - (distvalue1 + distvalue2)) < 1.e-6);
	assert(abs(data[1] - (distvalue2 * (distvalue3 - distvalue4) * // product
	                      distvalue2 * (distvalue5 - distvalue6))) < 1.e-6);

	// try to create a composite model / simulator where one of the submodels
	// does not implement one of the observables
	molstat::SimulateModelFactory subfactory_bad
		{ molstat::SimulateModelFactory::makeFactory<FailedSubModel>() };
	subfactory_bad.setDistribution("eps",
		make_shared<molstat::ConstantDistribution>(distvalue4));

	// try to make a simulator from the submodel... not allowed
	try
	{
		molstat::Simulator sim_bad{ subfactory_bad.getModel() };
		assert(false);
	}
	catch(const molstat::FullModelRequired &e)
	{
		// should be here
	}

	molstat::SimulateModelFactory cfactory_bad
		{ molstat::SimulateModelFactory::makeFactory<CompositeTestModelAdd>() };
	cfactory_bad.setDistribution("ef",
		make_shared<molstat::ConstantDistribution>(distvalue1));
	cfactory_bad.setDistribution("v",
		make_shared<molstat::ConstantDistribution>(distvalue2));
	cfactory_bad.addSubmodel(subfactory_bad.getModel());

	{
		molstat::Simulator sim{ cfactory_bad.getModel() };
		sim.setObservable(0, type_index{ typeid(BasicObs4) });
		data = sim.simulate(engine);
		assert(abs(data[0] - (distvalue1 + distvalue2)) < 1.e-6);

		// setting BasicObs1 should fail because the submodel doesn't support it
		try
		{
			sim.setObservable(1, type_index{ typeid(BasicObs1) });
			assert(false);
		}
		catch(const molstat::IncompatibleObservable &e)
		{
			// should be here
		}
	}

	return 0;
}
