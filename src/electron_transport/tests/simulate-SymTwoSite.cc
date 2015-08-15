/* This file is a part of MolStat, which is distributed under the Creative
   Commons Attribution-NonCommercial 4.0 International Public License.

   (c) 2014 Northwestern University. */

/**
 * \file tests/simulate-SymTwoSite.cc
 * \brief Test suite for the symmetric-coupling, two-site tight-binding
 *     channel.
 *
 * \test Test suite for the symmetric-coupling, two-site tight-binding
 *     channel.
 *
 * \author Matthew G.\ Reuter
 * \date October 2014
 */

#include <cassert>
#include <valarray>

#include <electron_transport/simulator_models/sym_two_site_channel.h>

using namespace std;

/// Shortcut for the type of channel used in this test.
using ChannelType = molstat::transport::SymTwoSiteChannel;

/**
 * \brief Main function for testing the symmetric-coupling, two-site model.
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
			.setDistribution("gamma", nullptr)
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
	auto SeebeckS = junction->getObservableFunction(
		type_index{ typeid(molstat::transport::SeebeckCoefficient) } );

	valarray<double> params(junction->get_num_parameters());

	// check known values for several parameter sets
	params[ChannelType::Index_EF] = 0.;
	params[ChannelType::Index_V] = 1.;
	params[ChannelType::Index_epsilon] = -4.;
	params[ChannelType::Index_gamma] = 0.8;
	params[ChannelType::Index_beta] = -3.;
	assert(abs(0.101007 - ZeroBiasG(params)) < thresh);
	assert(abs(0.127042 - ECurrent(params)) < thresh);
	assert(abs(0.127042 - StaticG(params)) < thresh);
	assert(abs(0.186815 - DiffG(params)) < thresh);
	assert(abs(params[ChannelType::Index_V] - AppBias(params)) < thresh);
	assert(abs(2.0089223 - SeebeckS(params)) < thresh);

	params[ChannelType::Index_EF] = 1.;
	params[ChannelType::Index_V] = -0.4;
	params[ChannelType::Index_epsilon] = -3.;
	params[ChannelType::Index_gamma] = 0.4;
	params[ChannelType::Index_beta] = -0.8;
	assert(abs(0.000431590 - ZeroBiasG(params)) < thresh);
	assert(abs(-0.000174208 - ECurrent(params)) < thresh);
	assert(abs(0.000435520 - StaticG(params)) < thresh);
	assert(abs(0.000443426 - DiffG(params)) < thresh);
	assert(abs(params[ChannelType::Index_V] - AppBias(params)) < thresh);
	assert(abs(1.0385126 - SeebeckS(params)) < thresh);

	params[ChannelType::Index_EF] = 3.;
	params[ChannelType::Index_V] = 1.4;
	params[ChannelType::Index_epsilon] = 1.1;
	params[ChannelType::Index_gamma] = 0.67;
	params[ChannelType::Index_beta] = -1.6;
	assert(abs(0.459683 - ZeroBiasG(params)) < thresh);
	assert(abs(0.673107 - ECurrent(params)) < thresh);
	assert(abs(0.480791 - StaticG(params)) < thresh);
	assert(abs(0.294527 - DiffG(params)) < thresh);
	assert(abs(params[ChannelType::Index_V] - AppBias(params)) < thresh);
	assert(abs(3.5332333 - SeebeckS(params)) < thresh);

	return 0;
}
