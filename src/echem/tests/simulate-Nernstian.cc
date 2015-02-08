/* This file is a part of MolStat, which is distributed under the Creative
   Commons Attribution-NonCommercial 4.0 International Public License.
   MolStat (c) 2014, Northwestern University. */
/**
 * \file tests/simulate-Nernstian.cc
 * \brief Test suite for Nernstian model of single molecule electrochemistry.
 *
 * \test Test suite for Nernstian model of single molecule electrochemistry.
 *
 * \author Bo Fu
 * \date January 2015
 */

#include <cassert>
#include <valarray>

#include <echem/simulator_models/nernstian.h>

using namespace std;

/// Shortcut for the type of model used in this test.
using ModelType = molstat::echem::NernstianReaction;

/**
 * \brief Main function for testing the Nernstian reaction model.
 *
 * \param[in] argc The number of command-line arguments.
 * \param[in] argv The command-line arguments.
 * \return Exit status: 0 if the code passes the test, non-zero otherwise.
 */
int main(int argc, char **argv)
{
	const double thresh = 1.0e-4;

	// use the factory to create a model
	shared_ptr<ModelType> nernstian = dynamic_pointer_cast<ModelType>(
			molstat::SimulateModelFactory::makeFactory<ModelType>()
			.setDistribution("af", nullptr)
			.setDistribution("ab", nullptr)
			.setDistribution("eref", nullptr)
			.getModel()
		);

	// get the observable functions
	auto FPotential = nernstian->getObservableFunction(
		type_index{ typeid(molstat::echem::ForwardETPotential) } );
	auto BPotential = nernstian->getObservableFunction(
		type_index{ typeid(molstat::echem::BackwardETPotential) } );

	valarray<double> params(3);

	// check known values for several parameter sets
	params[ModelType::Index_Af] = 5.0e3;
	params[ModelType::Index_Ab] = 5.0e3;
	params[ModelType::Index_Eref] = 1.1;
	assert(abs(1.1 - FPotential(params)) / 1.1 < thresh);
	assert(abs(1.1 - BPotential(params)) / 1.1 < thresh);

	params[ModelType::Index_Af] = 8.0e5;
	params[ModelType::Index_Ab] = 5.0e3;
	params[ModelType::Index_Eref] = 1.1;
	assert(abs(6.17517 - FPotential(params)) / 6.17517 < thresh);
	assert(abs(6.17517 - BPotential(params)) / 6.17517 < thresh);

	params[ModelType::Index_Af] = 8.0e5;
	params[ModelType::Index_Ab] = 9.0e8;
	params[ModelType::Index_Eref] = 0.4;
	assert(abs(-6.62554 - FPotential(params)) / 6.62554 < thresh);
	assert(abs(-6.62554 - BPotential(params)) / 6.62554 < thresh);

	return 0;
}
