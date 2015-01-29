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
	const double thresh = 1.0e-5;

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
	params[ModelType::Index_Af] = 5.0e3;
	params[ModelType::Index_Ab] = 5.0e3;
	params[ModelType::Index_Eref] = 1.1; 
	params[ModelType::Index_E0] = 0.;
	params[ModelType::Index_v] = 1.0;
	params[ModelType::Index_tlim] = 2.5;
	assert(abs(0.121622 - nonnernstian->E_applied(0.121622, params)) < thresh);
	assert(abs(2.32105 - nonnernstian->E_applied(2.678949, params)) < thresh);
	assert(abs(4.951881e3 - nonnernstian->kf(0.5, params)) < thresh * 1.0e4);
	assert(abs(4.123105e-7 - nonnernstian->kb(0.5, params)) < thresh * 1.0e-5);
	assert(abs(4.356003e-4 - nonnernstian->kf(1.5, params)) < thresh * 1.0e-3);
	assert(abs(2.284469e3 - nonnernstian->kb(1.5, params)) < thresh * 1.0e4);
	assert(abs(4.356003e-4 - nonnernstian->kf(3.5, params)) < thresh * 1.0e-3);
	assert(abs(2.284469e3 - nonnernstian->kb(3.5, params)) < thresh * 1.0e4);
	assert(abs(1.134077 - FPotential(params)) < thresh);
	assert(abs(1.065924 - BPotential(params)) < thresh);


	params[ModelType::Index_lambda] = 0.725;
	params[ModelType::Index_Af] = 5.0e6;
	params[ModelType::Index_Ab] = 5.0e6;
	params[ModelType::Index_Eref] = 1.2; 
	params[ModelType::Index_E0] = 0.;
	params[ModelType::Index_v] = 1.0e+3;
	params[ModelType::Index_tlim] = 2.5e-3;
	assert(abs(1.266938 - FPotential(params)) < thresh);
	assert(abs(1.133060 - BPotential(params)) < thresh);

	params[ModelType::Index_lambda] = 0.525;
	params[ModelType::Index_Af] = 5.0e9;
	params[ModelType::Index_Ab] = 5.0e9;
	params[ModelType::Index_Eref] = 0.9; 
	params[ModelType::Index_E0] = 0.;
	params[ModelType::Index_v] = 1.0e+6;
	params[ModelType::Index_tlim] = 2.5e-6;
	assert(abs(0.9150207 - FPotential(params)) < thresh);
	assert(abs(0.8849815 - BPotential(params)) < thresh);

	params[ModelType::Index_lambda] = 0.825;
	params[ModelType::Index_Af] = 5.0e12;
	params[ModelType::Index_Ab] = 5.0e12;
	params[ModelType::Index_Eref] = 1.0; 
	params[ModelType::Index_E0] = 0.;
	params[ModelType::Index_v] = 1.0e+9;
	params[ModelType::Index_tlim] = 2.5e-9;
	assert(abs(1.112267 - FPotential(params)) < thresh);
	assert(abs(0.8877257 - BPotential(params)) < thresh);


	return 0;
}
