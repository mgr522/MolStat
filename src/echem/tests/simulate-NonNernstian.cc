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
	params[ModelType::Index_lambda] = 24.176082;
	params[ModelType::Index_Af] = 5.0e3;
	params[ModelType::Index_Ab] = 5.0e3;
	params[ModelType::Index_Eref] = 42.549904;//reduced energy unit 
	params[ModelType::Index_E0] = 0.;//reduced energy unit 
	params[ModelType::Index_v] = 38.681731;//reduced energy unit
	params[ModelType::Index_tlim] = 2.5;
	assert(abs(4.704549 - nonnernstian->E_applied(0.121622, params)) < thresh);
	assert(abs(89.78227 - nonnernstian->E_applied(2.678949, params)) < thresh);
	assert(abs(4.356003e-4 - nonnernstian->kf(1.5, params)) < thresh * 1.0e-3);
	assert(abs(2.284469e3 - nonnernstian->kb(1.5, params)) < thresh * 1.0e4);
	assert(abs(4.356003e-4 - nonnernstian->kf(3.5, params)) < thresh * 1.0e-3);
	assert(abs(2.284469e3 - nonnernstian->kb(3.5, params)) < thresh * 1.0e4);
	assert(abs( (43.868061 - FPotential(params)) / 43.868061) < thresh);
	assert(abs( (41.23170 - BPotential(params)) / 41.23170) < thresh);


	params[ModelType::Index_lambda] = 28.044255;//reduced energy unit
	params[ModelType::Index_Af] = 5.0e6;
	params[ModelType::Index_Ab] = 5.0e6;
	params[ModelType::Index_Eref] = 46.418077; //reduced energy unit
	params[ModelType::Index_E0] = 0.;//reduced energy unit
	params[ModelType::Index_v] = 38.681731e+3;//reduced energy unit
	params[ModelType::Index_tlim] = 2.5e-3;
	assert(abs( (49.0073684 - FPotential(params)) / 49.0073684) < thresh);
	assert(abs( (43.828604 - BPotential(params)) / 43.828604) < thresh);

	params[ModelType::Index_lambda] = 20.307909;//reduced energy unit
	params[ModelType::Index_Af] = 5.0e9;
	params[ModelType::Index_Ab] = 5.0e9;
	params[ModelType::Index_Eref] = 34.813558; //reduced energy unit
	params[ModelType::Index_E0] = 0.;//reduced energy unit
	params[ModelType::Index_v] = 38.681731e+6;//reduced energy unit
	params[ModelType::Index_tlim] = 2.5e-6;
	assert(abs( (35.3945937 - FPotential(params)) / 35.3945937) < thresh);
	assert(abs( (34.2325201 - BPotential(params)) / 34.2325201) < thresh);

	params[ModelType::Index_lambda] = 31.912428;//reduced energy unit
	params[ModelType::Index_Af] = 5.0e12;
	params[ModelType::Index_Ab] = 5.0e12;
	params[ModelType::Index_Eref] = 38.681731; //reduced energy unit
	params[ModelType::Index_E0] = 0.;//reduced energy unit
	params[ModelType::Index_v] = 38.681731e+9;//reduced energy unit
	params[ModelType::Index_tlim] = 2.5e-9;
	assert(abs( (43.0244473 - FPotential(params)) / 43.0244473) < thresh);
	assert(abs( (34.3387563 - BPotential(params))/ 34.3387563)  < thresh);

}
