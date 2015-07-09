/* This file is a part of MolStat, which is distributed under the Creative
   Commons Attribution-NonCommercial 4.0 International Public License.

   (c) 2014 Northwestern University. */

/**
 * \file composite_simulate_model.cc
 * \brief Implements the molstat::CompositeSimulateModel class for simulating
 *    histograms.
 *
 * \author Matthew G.\ Reuter
 * \date October 2014
 */

#include "simulate_model.h"
#include "simulator_exceptions.h"

namespace molstat {

std::list<std::shared_ptr<const SimulateModel>>
	CompositeSimulateModel::getSubmodels() const
{
	std::list<std::shared_ptr<const SimulateModel>> ret;

	// go through the submodels and strip out the routing information
	for(const auto submodel : submodels)
		ret.emplace_back(submodel.first);

	return ret;
}

CompositeSimulateModel::SubmodelParameters
CompositeSimulateModel::routeSubmodelParameters(
	const std::valarray<double> &cparams) const
{
	SubmodelParameters ret;

	// go through all of the submodels
	// submodel.first is the pointer to the submodel
	// submodel.second is the list of indices
	for(const auto submodel : submodels)
		ret.emplace_back(std::make_pair(
			submodel.first,
			cparams[submodel.second]));

	return ret;
}

std::size_t CompositeSimulateModel::get_num_composite_parameters() const
{
	return get_names().size();
}

std::size_t CompositeSimulateModel::get_num_parameters() const
{
	// first get the number of parameters required directly
	std::size_t ret{ get_num_composite_parameters() };

	// add in the parameters for each submodel
	for(const auto submodel : submodels)
		ret += submodel.first->get_num_parameters();

	return ret;
}

std::valarray<double> CompositeSimulateModel::generateParameters(
	Engine &engine) const
{
	std::size_t tally = get_num_composite_parameters();
	std::valarray<double> ret(get_num_parameters());

	// simulate the parameters for the composite model
	for(std::size_t j = 0; j < tally; ++j)
	{
		ret[j] = dists[j]->sample(engine);
	}

	// go through the submodels, having them simulate their respective parameters
	for(const auto submodel : submodels)
	{
		std::size_t submodel_length = submodel.first->get_num_parameters();

		ret[std::slice(tally, submodel_length, 1)]
			= submodel.first->generateParameters(engine);

		// move the tally index up for the next model
		tally += submodel_length;
	}

	return ret;
}

} // namespace molstat
