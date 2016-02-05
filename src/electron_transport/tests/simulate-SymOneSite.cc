/* This file is a part of MolStat, which is distributed under the Creative
   Commons Attribution-NonCommercial 4.0 International Public License.

   (c) 2014 Northwestern University. */

/**
 * \file tests/simulate-SymOneSite.cc
 * \brief Test suite for the symmetric-coupling, single-site tight-binding
 *     channel.
 *
 * \test Test suite for the symmetric-coupling, single-site tight-binding
 *     channel.
 *
 * \author Matthew G.\ Reuter
 * \date October 2014
 */

#include <cassert>
#include <valarray>

#include <electron_transport/simulator_models/sym_one_site_channel.h>

using namespace std;

/// Shortcut for the type of channel used in this test.
using ChannelType = molstat::transport::SymOneSiteChannel;

/**
 * \brief Main function for testing the symmetric-coupling, one-site model.
 *
 * \param[in] argc The number of command-line arguments.
 * \param[in] argv The command-line arguments.
 * \return Exit status: 0 if the code passes the test, non-zero otherwise.
 */
int main(int argc, char **argv)
{
	const double thresh{ 1.e-5 };

	// use the factory to create a channel
	shared_ptr<ChannelType> channel = dynamic_pointer_cast<ChannelType>(
			molstat::SimulateModelFactory::makeFactory<ChannelType>()
			// not going to need the distributions, but the factory framework
			// requires them
			.setDistribution("epsilon", nullptr)
			.setDistribution("gamma", nullptr)
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
	auto ZeroBiasS = junction->getObservableFunction(
		type_index{ typeid(molstat::transport::ZeroBiasThermopower) } );

	valarray<double> params(junction->get_num_parameters());

	// check known values for several parameter sets
	params[ChannelType::Index_EF] = 0.;
	params[ChannelType::Index_V] = 1.;
	params[ChannelType::Index_epsilon] = -4.;
	params[ChannelType::Index_gamma] = 0.8;
	params[ChannelType::Index_a] = 0.;
	assert(abs(0.0384615 - ZeroBiasG(params)) < thresh);
	assert(abs(3.02309 - ECurrent(params)) < thresh);
	assert(abs(0.0390172 - StaticG(params)) < thresh);
	assert(abs(0.0401438 - DiffG(params)) < thresh);
	assert(abs(params[ChannelType::Index_V] - AppBias(params)) < thresh);
	assert(abs(0.480769 - ZeroBiasS(params)) < thresh);

	params[ChannelType::Index_EF] = 1.;
	params[ChannelType::Index_V] = -0.4;
	params[ChannelType::Index_epsilon] = -9.;
	params[ChannelType::Index_gamma] = 0.4;
	params[ChannelType::Index_a] = 0.;
	assert(abs(0.00159744 - ZeroBiasG(params)) < thresh);
	assert(abs(-0.0495283 - ECurrent(params)) < thresh);
	assert(abs(0.00159808 - StaticG(params)) < thresh);
	assert(abs(0.00159936 - DiffG(params)) < thresh);
	assert(abs(params[ChannelType::Index_V] - AppBias(params)) < thresh);
	assert(abs(0.199681 - ZeroBiasS(params)) < thresh);

	params[ChannelType::Index_EF] = 3.;
	params[ChannelType::Index_V] = 1.4;
	params[ChannelType::Index_epsilon] = -17.;
	params[ChannelType::Index_gamma] = 0.67;
	params[ChannelType::Index_a] = 0.;
	assert(abs(0.00112099 - ZeroBiasG(params)) < thresh);
	assert(abs(0.121746 - ECurrent(params)) < thresh);
	assert(abs(0.00112236 - StaticG(params)) < thresh);
	assert(abs(0.00112511 - DiffG(params)) < thresh);
	assert(abs(params[ChannelType::Index_V] - AppBias(params)) < thresh);
	assert(abs(0.0998879 - ZeroBiasS(params)) < thresh);


	params[ChannelType::Index_EF] = 0.;
	params[ChannelType::Index_V] = 1.;
	params[ChannelType::Index_epsilon] = -4.;
	params[ChannelType::Index_gamma] = 0.8;
	params[ChannelType::Index_a] = -0.1;
	assert(abs(0.0384615 - ZeroBiasG(params)) < thresh);
	assert(abs(2.88093 - ECurrent(params)) < thresh);
	assert(abs(0.0371825 - StaticG(params)) < thresh);
	assert(abs(0.0364382 - DiffG(params)) < thresh);
	assert(abs(params[ChannelType::Index_V] - AppBias(params)) < thresh);
	assert(abs(0.480769 - ZeroBiasS(params)) < thresh);


	params[ChannelType::Index_EF] = 0.;
	params[ChannelType::Index_V] = 1.;
	params[ChannelType::Index_epsilon] = 4.;
	params[ChannelType::Index_gamma] = 0.8;
	params[ChannelType::Index_a] = -0.1;
	assert(abs(0.0384615 - ZeroBiasG(params)) < thresh);
	assert(abs(3.17592 - ECurrent(params)) < thresh);
	assert(abs(0.0409897 - StaticG(params)) < thresh);
	assert(abs(0.0442754 - DiffG(params)) < thresh);
	assert(abs(params[ChannelType::Index_V] - AppBias(params)) < thresh);
	assert(abs(-0.480769 - ZeroBiasS(params)) < thresh);


	params[ChannelType::Index_EF] = 1.;
	params[ChannelType::Index_V] = -0.4;
	params[ChannelType::Index_epsilon] = -9.;
	params[ChannelType::Index_gamma] = 0.4;
	params[ChannelType::Index_a] = 1.;
	assert(abs(0.00159744 - ZeroBiasG(params)) < thresh);
	assert(abs(-0.0457959 - ECurrent(params)) < thresh);
	assert(abs(0.00147765 - StaticG(params)) < thresh);
	assert(abs(0.00136520 - DiffG(params)) < thresh);
	assert(abs(params[ChannelType::Index_V] - AppBias(params)) < thresh);
	assert(abs(0.199681 - ZeroBiasS(params)) < thresh);

	params[ChannelType::Index_EF] = 3.;
	params[ChannelType::Index_V] = 1.4;
	params[ChannelType::Index_epsilon] = -17.;
	params[ChannelType::Index_gamma] = 0.67;
	params[ChannelType::Index_a] = 0.24;
	assert(abs(0.00112099 - ZeroBiasG(params)) < thresh);
	assert(abs(0.125943 - ECurrent(params)) < thresh);
	assert(abs(0.00116105 - StaticG(params)) < thresh);
	assert(abs(0.00120367 - DiffG(params)) < thresh);
	assert(abs(params[ChannelType::Index_V] - AppBias(params)) < thresh);
	assert(abs(0.0998879 - ZeroBiasS(params)) < thresh);

	params[ChannelType::Index_EF] = 3.;
	params[ChannelType::Index_V] = 1.4;
	params[ChannelType::Index_epsilon] = -17.;
	params[ChannelType::Index_gamma] = 0.67;
	params[ChannelType::Index_a] = -0.05;
	assert(abs(0.00112099 - ZeroBiasG(params)) < thresh);
	assert(abs(0.120899 - ECurrent(params)) < thresh);
	assert(abs(0.00111445 - StaticG(params)) < thresh);
	assert(abs(0.00110948 - DiffG(params)) < thresh);
	assert(abs(params[ChannelType::Index_V] - AppBias(params)) < thresh);
	assert(abs(0.0998879 - ZeroBiasS(params)) < thresh);

	return 0;
}
