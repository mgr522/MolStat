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
	params[ModelType::Index_lambda] = 24.176082;//0.625eV
	params[ModelType::Index_Af] = 5.0e3;
	params[ModelType::Index_Ab] = 5.0e3;
	params[ModelType::Index_Eref] = 42.549904;//1.1V
	params[ModelType::Index_E0] = 0.;
	params[ModelType::Index_v] = 38.681731;//1.0V/s
	params[ModelType::Index_tlim] = 2.5;
	assert(abs( (58.02260 - nonnernstian->E_applied(1.5, params)) / 58.02260) < thresh);
	assert(abs( (4.35592e-4 - nonnernstian->kf(1.5, params)) / 4.35592e-4) < thresh);
	assert(abs( (2.284469e3 - nonnernstian->kb(1.5, params)) / 2.284469e3) < thresh);
	assert(abs( (54.15442 - nonnernstian->E_applied(3.6, params)) / 54.15442) < thresh);
	assert(abs( (8.90067e-3 - nonnernstian->kf(3.6, params)) / 8.90067e-3) < thresh);
	assert(abs( (9.75441e2 - nonnernstian->kb(3.6, params)) / 9.75441e2) < thresh);
	assert(abs( (43.868061 - FPotential(params)) / 43.868061) < thresh);
	assert(abs( (41.23178 - BPotential(params)) / 41.23178) < thresh);

	params[ModelType::Index_lambda] = 28.044255;//0.725eV
	params[ModelType::Index_Af] = 5.0e6;
	params[ModelType::Index_Ab] = 5.0e6;
	params[ModelType::Index_Eref] = 46.418077; //1.2V
	params[ModelType::Index_E0] = 0.;
	params[ModelType::Index_v] = 38.681731e+3;//1.0e+3
	params[ModelType::Index_tlim] = 2.5e-3;
	assert(abs( (49.0073684 - FPotential(params)) / 49.0073684) < thresh);
	assert(abs( (43.828604 - BPotential(params)) / 43.828604) < thresh);

	params[ModelType::Index_lambda] = 20.307909;//0.525
	params[ModelType::Index_Af] = 5.0e9;
	params[ModelType::Index_Ab] = 5.0e9;
	params[ModelType::Index_Eref] = 34.813558; //0.9
	params[ModelType::Index_E0] = 0.;
	params[ModelType::Index_v] = 38.681731e+6;//1.0e+6
	params[ModelType::Index_tlim] = 2.5e-6;
	assert(abs( (35.3945937 - FPotential(params)) / 35.3945937) < thresh);
	assert(abs( (34.2325201 - BPotential(params)) / 34.2325201) < thresh);

	params[ModelType::Index_lambda] = 31.912428;//0.825eV
	params[ModelType::Index_Af] = 5.0e12;
	params[ModelType::Index_Ab] = 5.0e12;
	params[ModelType::Index_Eref] = 38.681731; //1.0V
	params[ModelType::Index_E0] = 0.;
	params[ModelType::Index_v] = 38.681731e+9;//1.0e+9
	params[ModelType::Index_tlim] = 2.5e-9;
	assert(abs( (43.0244473 - FPotential(params)) / 43.0244473) < thresh);
	assert(abs( (34.3387563 - BPotential(params))/ 34.3387563)  < thresh);

}
