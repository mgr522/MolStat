/* This file is a part of MolStat, which is distributed under the Creative
   Commons Attribution-NonCommercial 4.0 International Public License.

   (c) 2014 Northwestern University. */
/**
 * \file simulate_model_factory.cc
 * \brief Implements the molstat::SimulateModelFactory class for constructing
 *    molstat::SimulateModel objects.
 *
 * \author Matthew G.\ Reuter
 * \date October 2014
 */

#include "simulate_model.h"
#include "simulator_exceptions.h"

namespace molstat {

SimulateModelFactory &SimulateModelFactory::setDistribution(std::string name,
	std::shared_ptr<const RandomDistribution> dist,
	bool *used_dist /* = nullptr */)
{
	// set the distribution in the model
	const std::size_t length = model->get_num_parameters();
	const std::string lower_name{ to_lower(name) };

	if(used_dist != nullptr)
		*used_dist = false;

	// find the position of the distribution
	for(std::size_t pos = 0; pos < length; ++pos)
	{
		if(lower_name == to_lower(model_names[pos]))
		{
			// found! set the distribution and signal that the distribution is used
			model->dists[pos] = dist;
			if(used_dist != nullptr)
				*used_dist = true;

			// check off the name
			auto spot = remaining_names.find(lower_name);
			if(spot != remaining_names.end())
				remaining_names.erase(spot);
		}
	}

	return *this;
}

SimulateModelFactory &SimulateModelFactory::addSubmodel(
	std::shared_ptr<SimulateModel> submodel_add)
{
	// make sure this is a composite model
	if(comp_model == nullptr)
		throw NotCompositeSimulateModel();

	// make sure the submodel is of the correct type
	if(submodel_add->getModelType() != comp_model->getSubmodelType())
		throw IncompatibleSubmodel();

	// add the model
	comp_model->submodels.emplace_back(submodel_add);

	return *this;
}

std::shared_ptr<SimulateModel> SimulateModelFactory::getModel()
{
	// make sure all distributions have been specified
	if(remaining_names.size() > 0)
		throw MissingDistribution(*(remaining_names.cbegin()));

	// if this is a composite model, make sure at least one submodel has been
	// specified
	if(comp_model != nullptr && comp_model->submodels.size() == 0)
		throw NoSubmodels();

	return model;
}

} // namespace molstat
