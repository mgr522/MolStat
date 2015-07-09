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

const std::list<std::shared_ptr<SimulateModel>> &
	CompositeSimulateModel::getSubmodels() const
{
	return submodels;
CompositeSimulateModel::SubmodelParameters
CompositeSimulateModel::routeSubmodelParameters(
	const std::valarray<double> &cparams) const
{
	SubmodelParameters ret;

	// get the number of model parameters explicitly required by the
	// composite model
	const size_t cnparam{ get_num_composite_parameters() };
	size_t tally { cnparam };

	// go through all of the submodels
	for(const auto submodel : submodels)
	{
		const size_t sub_nparam{ submodel->get_num_parameters() };

		// create a list of the parameters indices that should be passed
		// to this submodel
		std::valarray<size_t> indices(cnparam + sub_nparam);
		// first add in the indices required by the composite model
		for(size_t j = 0; j < cnparam; ++j)
			indices[j] = j;
		// now add in the indices for the specific submodel
		for(size_t j = tally; j < tally + sub_nparam; ++j)
			indices[j - tally + cnparam] = j;

		ret.emplace_back(make_pair(
			submodel,
			cparams[indices]));

		// now add the submodel's parameters to the tally for the next offset
		tally += sub_nparam;
	}

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
		ret += submodel->get_num_parameters();

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

	// go through the submodels, having them simulate parameters
	for(const auto submodel : submodels)
	{
		std::size_t submodel_length = submodel->get_num_parameters();

		ret[std::slice(tally, submodel_length, 1)]
			= submodel->generateParameters(engine);

		// move the tally index up for the next model
		tally += submodel_length;
	}

	return ret;
}

} // namespace molstat
