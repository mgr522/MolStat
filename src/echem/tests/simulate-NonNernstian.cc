/* This file is a part of MolStat, which is distributed under the Creative
   Commons Attribution-NonCommercial 4.0 International Public License.
   MolStat (c) 2014, Northwestern University. */
/**
 * \file tests/simulate-NonNernstian.cc
 * \brief Test suite for non-Nernstian model of single molecule
 *    electrochemistry.
 *
 * \test Test suite for non-Nernstian model of single molecule
 *    electrochemistry.
 *
 * \author Bo Fu
 * \date January 2015
 */

#include <cassert>
#include <valarray>
#include <iostream>

#include <echem/simulator_models/non-nernstian.h>

using namespace std;

/// Shortcut for the type of channel used in this test.
using ModelType = molstat::echem::NonNernstianReaction;

/**
 * \brief Main function for testing the non-Nernstian reaction model.
 *
 * \param[in] argc The number of command-line arguments.
 * \param[in] argv The command-line arguments.
 * \return Exit status: 0 if the code passes the test, non-zero otherwise.
 */
int main(int argc, char **argv)
{
	const double thresh = 1.0e-4;

	// use the factory to create a junction
	shared_ptr<ModelType> nonnernstian = dynamic_pointer_cast<ModelType>(
			molstat::SimulateModelFactory::makeFactory<ModelType>()
			.setDistribution("lambda", nullptr)
			.setDistribution("af", nullptr)
			.setDistribution("ab", nullptr)
			.setDistribution("eref", nullptr)
			.setDistribution("e0", nullptr)
			.setDistribution("v", nullptr)
			.setDistribution("tlim", nullptr)
			.getModel()
		);

	// get the observable functions
	auto FPotential = nonnernstian->getObservableFunction(
		type_index{ typeid(molstat::echem::ForwardETPotential) } );
	auto BPotential = nonnernstian->getObservableFunction(
		type_index{ typeid(molstat::echem::BackwardETPotential) } );

	valarray<double> params(7);

	// check known values for several parameter sets
	params[ModelType::Index_lambda] = 0.625;
	params[ModelType::Index_Af] = 5.e3;
	params[ModelType::Index_Ab] = 5.e3;
	params[ModelType::Index_Eref] = 1.1;
	params[ModelType::Index_E0] = 0.;
	params[ModelType::Index_v] = 1.;
	params[ModelType::Index_tlim] = 2.5;
	assert(abs((1.10012 - FPotential(params)) / 1.10012) < thresh);
	assert(abs((1.09988 - BPotential(params)) / 1.09988) < thresh);

	params[ModelType::Index_lambda] = 0.725;
	params[ModelType::Index_Af] = 5.e6;
	params[ModelType::Index_Ab] = 5.e6;
	params[ModelType::Index_Eref] = 1.2;
	params[ModelType::Index_E0] = 0.;
	params[ModelType::Index_v] = 1.e3;
	params[ModelType::Index_tlim] = 2.5e-3;
	assert(abs((1.20010 - FPotential(params)) / 1.20010) < thresh);
	assert(abs((1.19987 - BPotential(params)) / 1.19987) < thresh);

	params[ModelType::Index_lambda] = 0.525;
	params[ModelType::Index_Af] = 5.e9;
	params[ModelType::Index_Ab] = 5.e9;
	params[ModelType::Index_Eref] = 0.9;
	params[ModelType::Index_E0] = 0.;
	params[ModelType::Index_v] = 1.e6;
	params[ModelType::Index_tlim] = 2.5e-6;
	assert(abs((0.900104 - FPotential(params)) / 0.900104) < thresh);
	assert(abs((0.899887 - BPotential(params)) / 0.899887) < thresh);

	params[ModelType::Index_lambda] = 0.825;
	params[ModelType::Index_Af] = 5.e12;
	params[ModelType::Index_Ab] = 5.e12;
	params[ModelType::Index_Eref] = 1.;
	params[ModelType::Index_E0] = 0.;
	params[ModelType::Index_v] = 1.e9;
	params[ModelType::Index_tlim] = 2.5e-9;
	assert(abs((1.00012 - FPotential(params)) / 1.00012) < thresh);
	assert(abs((0.999889 - BPotential(params))/ 0.999889)  < thresh);

	// this test yields a system where the desired potentials are not found...
	// make sure the exception is thrown
	params[ModelType::Index_lambda] = 0.825;
	params[ModelType::Index_Af] = 5.e6;
	params[ModelType::Index_Ab] = 5.e4;
	params[ModelType::Index_Eref] = 1.;
	params[ModelType::Index_E0] = -1.;
	params[ModelType::Index_v] = 1.;
	params[ModelType::Index_tlim] = 2.5;
	try
	{
		FPotential(params);
		assert(false); // should have thrown
	}
	catch(const molstat::NoObservableProduced &e)
	{
		// should be here
	}
	try
	{
		BPotential(params);
		assert(false); // should have thrown
	}
	catch(const molstat::NoObservableProduced &e)
	{
		// should be here
	}

	return 0;
}
