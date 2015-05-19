/* This file is a part of MolStat, which is distributed under the Creative
   Commons Attribution-NonCommercial 4.0 International Public License.

   (c) 2015 Stony Brook University. */

/**
 * \file tests/simulate-RectBarrier.cc
 * \brief Test suite for the rectangular barrier channel model.
 *
 * \test Test suite for the rectangular barrier channel model.
 *
 * \author Matthew G.\ Reuter
 * \date January 2015
 */

#include <cassert>
#include <valarray>
#include <iostream>

#include <electron_transport/simulator_models/rectangular_barrier.h>

using namespace std;

/// Shortcut for the type of channel used in this test.
using ChannelType = molstat::transport::RectangularBarrier;

/**
 * \brief Main function for testing the symmetric-coupling, one-site model.
 *
 * \param[in] argc The number of command-line arguments.
 * \param[in] argv The command-line arguments.
 * \return Exit status: 0 if the code passes the test, non-zero otherwise.
 */
int main(int argc, char **argv)
{
	const double thresh{ 1.e-6 };

	// use the factory to create a channel
	shared_ptr<ChannelType> channel = dynamic_pointer_cast<ChannelType>(
			molstat::SimulateModelFactory::makeFactory<ChannelType>()
			// not going to need the distributions, but the factory framework
			// requires them
			.setDistribution("height", nullptr)
			.setDistribution("width", nullptr)
			.getModel()
		);
	assert(channel != nullptr);

	// use the factory to create a junction
	shared_ptr<molstat::SimulateModel> junction =
		molstat::SimulateModelFactory::makeFactory
			<molstat::transport::TransportJunction>()
		.setDistribution("ef", nullptr)
		.setDistribution("v", nullptr)
		.addSubmodel(channel)
		.getModel();

	// get the observable functions
	auto AppBias = junction->getObservableFunction(
		type_index{ typeid(molstat::transport::AppliedBias) } );
	auto ZeroBiasG = junction->getObservableFunction(
		type_index{ typeid(molstat::transport::ZeroBiasConductance) } );
	
	valarray<double> params(junction->get_num_parameters());

	// check known values for several parameter sets
	params[ChannelType::Index_EF] = 0.2;
	params[ChannelType::Index_V] = 0.;
	params[ChannelType::Index_h] = 1.4;
	params[ChannelType::Index_w] = 1.;
	assert(abs(2.61472e-5 - ZeroBiasG(params)) < thresh);
	assert(abs(params[ChannelType::Index_V] - AppBias(params)) < thresh);

	params[ChannelType::Index_EF] = 0.3;
	params[ChannelType::Index_V] = 0.;
	params[ChannelType::Index_h] = 0.6;
	params[ChannelType::Index_w] = 0.5;
	assert(abs(0.214993 - ZeroBiasG(params)) < thresh);
	assert(abs(params[ChannelType::Index_V] - AppBias(params)) < thresh);

	return 0;
}
