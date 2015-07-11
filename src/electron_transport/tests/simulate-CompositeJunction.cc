/* This file is a part of MolStat, which is distributed under the Creative
   Commons Attribution-NonCommercial 4.0 International Public License.

   (c) 2014 Northwestern University. */

/**
 * \file tests/simulate-CompositeJunction.cc
 * \brief Test suite for a composite junction with two channels.
 *
 * \test Test suite for a composite junction with two channels.
 *
 * \author Matthew G.\ Reuter
 * \date October 2014
 */

#include <cassert>
#include <valarray>

#include <electron_transport/simulator_models/sym_one_site_channel.h>
#include <electron_transport/simulator_models/asym_one_site_channel.h>

using namespace std;

/**
 * \brief Test a composite junction with one symmetrically-coupled
 *    and one asymmetrically-coupled channel.
 *
 * \param[in] thresh The numerical threshold for assessing the tests.
 */
static void test_sym_asym(const double thresh) {
	// Alias for one type of channel used in this test.
	using ChannelType1 = molstat::transport::SymOneSiteChannel;
	// Alias for the other type of channel used in this test.
	using ChannelType2 = molstat::transport::AsymOneSiteChannel;

	// use the factory to create the channel
	shared_ptr<ChannelType1> channel1 = dynamic_pointer_cast<ChannelType1>(
			molstat::SimulateModelFactory::makeFactory<ChannelType1>()
			// not going to need the distributions, but the factory framework
			// requires them
			.setDistribution("epsilon", nullptr)
			.setDistribution("gamma", nullptr)
			.setDistribution("a", nullptr)
			.getModel()
		);
	assert(channel1 != nullptr);

	shared_ptr<ChannelType2> channel2 = dynamic_pointer_cast<ChannelType2>(
			molstat::SimulateModelFactory::makeFactory<ChannelType2>()
			.setDistribution("epsilon", nullptr)
			.setDistribution("gammal", nullptr)
			.setDistribution("gammar", nullptr)
			.setDistribution("a", nullptr)
			.getModel()
		);
	assert(channel2 != nullptr);

	// use the factory to create a junction
	shared_ptr<molstat::SimulateModel> junction =
		molstat::SimulateModelFactory::makeFactory
			<molstat::transport::TransportJunction>()
		.setDistribution("ef", nullptr)
		.setDistribution("v", nullptr)
		.addSubmodel(channel1)
		.addSubmodel(channel2)
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

	// subarrays for passing parameters directly to the submodels
	const valarray<size_t> index1 = {0, 1, 2, 3, 4};
	const valarray<size_t> index2 = {0, 1, 5, 6, 7, 8};

	// check known values for several parameter sets
	params[ChannelType1::Index_EF] = 0.;
	params[ChannelType1::Index_V] = 1.;
	params[ChannelType1::Index_epsilon] = 3.4;
	params[ChannelType1::Index_gamma] = 0.4;
	params[ChannelType1::Index_a] = 0.04;
	params[ChannelType2::Index_epsilon + channel1->get_num_parameters()] = -4.;
	params[ChannelType2::Index_gammaL + channel1->get_num_parameters()] = 0.8;
	params[ChannelType2::Index_gammaR + channel1->get_num_parameters()] = 1.;
	params[ChannelType2::Index_a + channel1->get_num_parameters()] = 0.;
	assert(abs(channel1->ECurrent(params[index1]) +
		channel2->ECurrent(params[index2]) - ECurrent(params)) < thresh);
	assert(abs(channel1->ZeroBiasG(params[index1]) +
		channel2->ZeroBiasG(params[index2]) - ZeroBiasG(params)) < thresh);
	assert(abs(channel1->StaticG(params[index1]) +
		channel2->StaticG(params[index2]) - StaticG(params)) < thresh);
	assert(abs(channel1->DiffG(params[index1]) + channel2->DiffG(params[index2])
		- DiffG(params)) < thresh);
	assert(abs(params[ChannelType1::Index_V] - AppBias(params)) < thresh);

	params[ChannelType1::Index_EF] = 0.;
	params[ChannelType1::Index_V] = 1.;
	params[ChannelType1::Index_epsilon] = -4.1;
	params[ChannelType1::Index_gamma] = 0.7;
	params[ChannelType1::Index_a] = -0.04;
	params[ChannelType2::Index_epsilon + channel1->get_num_parameters()] = -9.;
	params[ChannelType2::Index_gammaL + channel1->get_num_parameters()] = 0.4;
	params[ChannelType2::Index_gammaR + channel1->get_num_parameters()] = 0.2;
	params[ChannelType2::Index_a + channel1->get_num_parameters()] = 0.;
	assert(abs(channel1->ECurrent(params[index1]) +
		channel2->ECurrent(params[index2]) - ECurrent(params)) < thresh);
	assert(abs(channel1->ZeroBiasG(params[index1]) +
		channel2->ZeroBiasG(params[index2]) - ZeroBiasG(params)) < thresh);
	assert(abs(channel1->StaticG(params[index1]) +
		channel2->StaticG(params[index2]) - StaticG(params)) < thresh);
	assert(abs(channel1->DiffG(params[index1]) + channel2->DiffG(params[index2])
		- DiffG(params)) < thresh);
	assert(abs(params[ChannelType1::Index_V] - AppBias(params)) < thresh);

	params[ChannelType1::Index_EF] = 0.;
	params[ChannelType1::Index_V] = 1.;
	params[ChannelType1::Index_epsilon] = 0.6;
	params[ChannelType1::Index_gamma] = 0.3;
	params[ChannelType1::Index_a] = 0.;
	params[ChannelType2::Index_epsilon + channel1->get_num_parameters()] = -17.;
	params[ChannelType2::Index_gammaL + channel1->get_num_parameters()] = 0.67;
	params[ChannelType2::Index_gammaR + channel1->get_num_parameters()] = 1.98;
	params[ChannelType2::Index_a + channel1->get_num_parameters()] = 0.;
	assert(abs(channel1->ECurrent(params[index1]) +
		channel2->ECurrent(params[index2]) - ECurrent(params)) < thresh);
	assert(abs(channel1->ZeroBiasG(params[index1]) +
		channel2->ZeroBiasG(params[index2]) - ZeroBiasG(params)) < thresh);
	assert(abs(channel1->StaticG(params[index1]) +
		channel2->StaticG(params[index2]) - StaticG(params)) < thresh);
	assert(abs(channel1->DiffG(params[index1]) + channel2->DiffG(params[index2])
		- DiffG(params)) < thresh);
	assert(abs(params[ChannelType1::Index_V] - AppBias(params)) < thresh);

	params[ChannelType1::Index_EF] = 0.;
	params[ChannelType1::Index_V] = 1.;
	params[ChannelType1::Index_epsilon] = 3.4;
	params[ChannelType1::Index_gamma] = 3.;
	params[ChannelType1::Index_a] = -0.03;
	params[ChannelType2::Index_epsilon + channel1->get_num_parameters()] = -4.;
	params[ChannelType2::Index_gammaL + channel1->get_num_parameters()] = 0.8;
	params[ChannelType2::Index_gammaR + channel1->get_num_parameters()] = 1.;
	params[ChannelType2::Index_a + channel1->get_num_parameters()] = 0.1;
	assert(abs(channel1->ECurrent(params[index1]) +
		channel2->ECurrent(params[index2]) - ECurrent(params)) < thresh);
	assert(abs(channel1->ZeroBiasG(params[index1]) +
		channel2->ZeroBiasG(params[index2]) - ZeroBiasG(params)) < thresh);
	assert(abs(channel1->StaticG(params[index1]) +
		channel2->StaticG(params[index2]) - StaticG(params)) < thresh);
	assert(abs(channel1->DiffG(params[index1]) + channel2->DiffG(params[index2])
		- DiffG(params)) < thresh);
	assert(abs(params[ChannelType1::Index_V] - AppBias(params)) < thresh);

	// call to the Seebeck coefficient routines should throw because the
	// asymmetric model does not presently support the Seebeck coefficient.
	try {
		SeebeckS(params);
		assert(false);
	}
	catch(const molstat::IncompatibleObservable &e) {
		// should be here
	}
	catch(...) {
		assert(false);
	}
}

/**
 * \brief Test a composite junction with one symmetrically-coupled
 *    and one symmetrically-coupled channel.
 *
 * \param[in] thresh The numerical threshold for assessing the tests.
 */
static void test_sym_sym(const double thresh) {
	// Alias for one type of channel used in this test.
	using ChannelType1 = molstat::transport::SymOneSiteChannel;
	// Alias for the other type of channel used in this test.
	using ChannelType2 = molstat::transport::SymOneSiteChannel;

	// use the factory to create the channel
	shared_ptr<ChannelType1> channel1 = dynamic_pointer_cast<ChannelType1>(
			molstat::SimulateModelFactory::makeFactory<ChannelType1>()
			// not going to need the distributions, but the factory framework
			// requires them
			.setDistribution("epsilon", nullptr)
			.setDistribution("gamma", nullptr)
			.setDistribution("a", nullptr)
			.getModel()
		);
	assert(channel1 != nullptr);

	shared_ptr<ChannelType2> channel2 = dynamic_pointer_cast<ChannelType2>(
			molstat::SimulateModelFactory::makeFactory<ChannelType2>()
			.setDistribution("epsilon", nullptr)
			.setDistribution("gamma", nullptr)
			.setDistribution("a", nullptr)
			.getModel()
		);
	assert(channel2 != nullptr);

	// use the factory to create a junction
	shared_ptr<molstat::SimulateModel> junction =
		molstat::SimulateModelFactory::makeFactory
			<molstat::transport::TransportJunction>()
		.setDistribution("ef", nullptr)
		.setDistribution("v", nullptr)
		.addSubmodel(channel1)
		.addSubmodel(channel2)
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

	// subarrays for passing parameters directly to the submodels
	const valarray<size_t> index1 = {0, 1, 2, 3, 4};
	const valarray<size_t> index2 = {0, 1, 5, 6, 7};

	// check known values for several parameter sets
	params[ChannelType1::Index_EF] = 0.;
	params[ChannelType1::Index_V] = 1.;
	params[ChannelType1::Index_epsilon] = 3.4;
	params[ChannelType1::Index_gamma] = 0.4;
	params[ChannelType1::Index_a] = 0.04;
	params[ChannelType2::Index_epsilon + channel1->get_num_parameters()] = -4.;
	params[ChannelType2::Index_gamma + channel1->get_num_parameters()] = 0.8;
	params[ChannelType2::Index_a + channel1->get_num_parameters()] = 0.;
	assert(abs(channel1->ECurrent(params[index1]) +
		channel2->ECurrent(params[index2]) - ECurrent(params)) < thresh);
	assert(abs(channel1->ZeroBiasG(params[index1]) +
		channel2->ZeroBiasG(params[index2]) - ZeroBiasG(params)) < thresh);
	assert(abs(channel1->StaticG(params[index1]) +
		channel2->StaticG(params[index2]) - StaticG(params)) < thresh);
	assert(abs(channel1->DiffG(params[index1]) + channel2->DiffG(params[index2])
		- DiffG(params)) < thresh);
	assert(abs(params[ChannelType1::Index_V] - AppBias(params)) < thresh);
	assert(abs(
		(channel1->SeebeckS(params[index1])*channel1->ZeroBiasG(params[index1])
			+channel2->SeebeckS(params[index2])*channel2->ZeroBiasG(params[index2]))
			/ ZeroBiasG(params) - SeebeckS(params)
		) < thresh);
}

/**
 * \brief Main function for testing the junction model with multiple channels.
 *
 * \param[in] argc The number of command-line arguments.
 * \param[in] argv The command-line arguments.
 * \return Exit status: 0 if the code passes the test, non-zero otherwise.
 */
int main(int argc, char **argv)
{
	const double thresh{ 1.e-6 };

	test_sym_sym(thresh);
	test_sym_asym(thresh);
	
	return 0;
}
