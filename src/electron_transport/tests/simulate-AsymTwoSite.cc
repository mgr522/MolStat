/* This file is a part of MolStat, which is distributed under the Creative
   Commons Attribution-NonCommercial 4.0 International Public License.

   (c) 2014 Northwestern University. */

/**
 * \file tests/simulate-AsymTwoSite.cc
 * \brief Test suite for the asymmetric-coupling, two-site tight-binding
 *     channel.
 *
 * \test Test suite for the asymmetric-coupling, two-site tight-binding
 *     channel.
 *
 * \author Matthew G.\ Reuter
 * \date October 2014
 */

#include <cassert>
#include <valarray>

#include <electron_transport/simulator_models/asym_two_site_channel.h>

using namespace std;

/// Shortcut for the type of channel used in this test.
using ChannelType = molstat::transport::AsymTwoSiteChannel;

/**
 * \brief Main function for testing the asymmetric-coupling, two-site model.
 *
 * \param[in] argc The number of command-line arguments.
 * \param[in] argv The command-line arguments.
 * \return Exit status: 0 if the code passes the test, non-zero otherwise.
 */
int main(int argc, char **argv)
{
	const double thresh = 1.0e-5;

	// use the factory to create a channel
	shared_ptr<ChannelType> channel = dynamic_pointer_cast<ChannelType>(
			molstat::SimulateModelFactory::makeFactory<ChannelType>()
			// not going to need the distributions, but the factory framework
			// requires them
			.setDistribution("epsilon", nullptr)
			.setDistribution("gammal", nullptr)
			.setDistribution("gammar", nullptr)
			.setDistribution("beta", nullptr)
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
	params[ChannelType::Index_beta] = -3.;
	assert(abs(0.121622 - ZeroBiasG(params)) < thresh);
	assert(abs(11.6171 - ECurrent(params)) / 11.6171 < thresh);
	assert(abs(0.149936 - StaticG(params)) < thresh);
	assert(abs(0.213248 - DiffG(params)) < thresh);
	assert(abs(params[ChannelType::Index_V] - AppBias(params)) < thresh);

	params[ChannelType::Index_EF] = 1.;
	params[ChannelType::Index_V] = -0.4;
	params[ChannelType::Index_epsilon] = -3.;
	params[ChannelType::Index_gammaL] = 0.4;
	params[ChannelType::Index_gammaR] = 0.2;
	params[ChannelType::Index_beta] = -0.8;
	assert(abs(0.000216257 - ZeroBiasG(params)) < thresh);
	assert(abs(-0.0067635 - ECurrent(params)) < thresh);
	assert(abs(0.000218231 - StaticG(params)) < thresh);
	assert(abs(0.000222203 - DiffG(params)) < thresh);
	assert(abs(params[ChannelType::Index_V] - AppBias(params)) < thresh);

	params[ChannelType::Index_EF] = -1.;
	params[ChannelType::Index_V] = 1.4;
	params[ChannelType::Index_epsilon] = 5.;
	params[ChannelType::Index_gammaL] = 0.67;
	params[ChannelType::Index_gammaR] = 1.98;
	params[ChannelType::Index_beta] = -1.6;
	assert(abs(0.00292927 - ZeroBiasG(params)) < thresh);
	assert(abs(0.3345 - ECurrent(params)) < thresh);
	assert(abs(0.00308371 - StaticG(params)) < thresh);
	assert(abs(0.00340305 - DiffG(params)) < thresh);
	assert(abs(params[ChannelType::Index_V] - AppBias(params)) < thresh);

	return 0;
}
