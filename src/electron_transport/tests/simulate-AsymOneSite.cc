/* This file is a part of MolStat, which is distributed under the Creative
   Commons Attribution-NonCommercial 4.0 International Public License.

   (c) 2014 Northwestern University. */

/**
 * \file tests/simulate-AsymOneSite.cc
 * \brief Test suite for the asymmetric-coupling, single-site tight-binding
 *     channel.
 *
 * \test Test suite for the asymmetric-coupling, single-site tight-binding
 *     channel.
 *
 * \author Matthew G.\ Reuter
 * \date October 2014
 */

#include <cassert>
#include <valarray>

#include <electron_transport/simulator_models/asym_one_site_channel.h>

using namespace std;

/// Shortcut for the type of channel used in this test.
using ChannelType = molstat::transport::AsymOneSiteChannel;

/**
 * \brief Main function for testing the asymmetric-coupling, one-site model.
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
			.setDistribution("epsilon", nullptr)
			.setDistribution("gammal", nullptr)
			.setDistribution("gammar", nullptr)
			.setDistribution("a", nullptr)
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
	auto ECurrent = junction->getObservableFunction(
		type_index{ typeid(molstat::transport::ElectricCurrent) } );
	auto ZeroBiasG = junction->getObservableFunction(
		type_index{ typeid(molstat::transport::ZeroBiasConductance) } );
	auto StaticG = junction->getObservableFunction(
		type_index{ typeid(molstat::transport::StaticConductance) } );
	auto DiffG = junction->getObservableFunction(
		type_index{ typeid(molstat::transport::DifferentialConductance) } );

	valarray<double> params(junction->get_num_parameters());

	// check known values for several parameter sets
	params[ChannelType::Index_EF] = 0.;
	params[ChannelType::Index_V] = 1.;
	params[ChannelType::Index_epsilon] = -4.;
	params[ChannelType::Index_gammaL] = 0.8;
	params[ChannelType::Index_gammaR] = 1.;
	params[ChannelType::Index_a] = 0.;
	assert(abs(0.0475907 - ZeroBiasG(params)) < thresh);
	assert(abs(0.0482617 - ECurrent(params)) < thresh);
	assert(abs(0.0482617 - StaticG(params)) < thresh);
	assert(abs(0.0496212 - DiffG(params)) < thresh);
	assert(abs(params[ChannelType::Index_V] - AppBias(params)) < thresh);

	params[ChannelType::Index_EF] = 1.;
	params[ChannelType::Index_V] = -0.4;
	params[ChannelType::Index_epsilon] = -9.;
	params[ChannelType::Index_gammaL] = 0.4;
	params[ChannelType::Index_gammaR] = 0.2;
	params[ChannelType::Index_a] = 0.;
	assert(abs(0.000799281 - ZeroBiasG(params)) < thresh);
	assert(abs(-0.000319840 - ECurrent(params)) < thresh);
	assert(abs(0.000799600 - StaticG(params)) < thresh);
	assert(abs(0.000800238 - DiffG(params)) < thresh);
	assert(abs(params[ChannelType::Index_V] - AppBias(params)) < thresh);

	params[ChannelType::Index_EF] = 3.;
	params[ChannelType::Index_V] = 1.4;
	params[ChannelType::Index_epsilon] = -17.;
	params[ChannelType::Index_gammaL] = 0.67;
	params[ChannelType::Index_gammaR] = 1.98;
	params[ChannelType::Index_a] = 0.;
	assert(abs(0.00330201 - ZeroBiasG(params)) < thresh);
	assert(abs(0.00462842 - ECurrent(params)) < thresh);
	assert(abs(0.00330602 - StaticG(params)) < thresh);
	assert(abs(0.00331404 - DiffG(params)) < thresh);
	assert(abs(params[ChannelType::Index_V] - AppBias(params)) < thresh);

	params[ChannelType::Index_EF] = 0.;
	params[ChannelType::Index_V] = 1.;
	params[ChannelType::Index_epsilon] = -4.;
	params[ChannelType::Index_gammaL] = 0.8;
	params[ChannelType::Index_gammaR] = 1.;
	params[ChannelType::Index_a] = 0.1;
	assert(abs(0.0475907 - ZeroBiasG(params)) < thresh);
	assert(abs(0.0506743 - ECurrent(params)) < thresh);
	assert(abs(0.0506743 - StaticG(params)) < thresh);
	assert(abs(0.0546687 - DiffG(params)) < thresh);
	assert(abs(params[ChannelType::Index_V] - AppBias(params)) < thresh);

	params[ChannelType::Index_EF] = 1.;
	params[ChannelType::Index_V] = -0.4;
	params[ChannelType::Index_epsilon] = -9.;
	params[ChannelType::Index_gammaL] = 0.4;
	params[ChannelType::Index_gammaR] = 0.2;
	params[ChannelType::Index_a] = -0.03;
	assert(abs(0.000799281 - ZeroBiasG(params)) < thresh);
	assert(abs(-0.000320609 - ECurrent(params)) < thresh);
	assert(abs(0.000801521 - StaticG(params)) < thresh);
	assert(abs(0.000804088 - DiffG(params)) < thresh);
	assert(abs(params[ChannelType::Index_V] - AppBias(params)) < thresh);

	params[ChannelType::Index_EF] = 3.;
	params[ChannelType::Index_V] = 1.4;
	params[ChannelType::Index_epsilon] = -17.;
	params[ChannelType::Index_gammaL] = 0.67;
	params[ChannelType::Index_gammaR] = 1.98;
	params[ChannelType::Index_a] = 1.3;
	assert(abs(0.00330201 - ZeroBiasG(params)) < thresh);
	assert(abs(0.00559778 - ECurrent(params)) < thresh);
	assert(abs(0.00399841 - StaticG(params)) < thresh);
	assert(abs(0.00480763 - DiffG(params)) < thresh);
	assert(abs(params[ChannelType::Index_V] - AppBias(params)) < thresh);

	return 0;
}
