/* This file is a part of MolStat, which is distributed under the Creative
   Commons Attribution-NonCommercial 4.0 International Public License.

   (c) 2014 Northwestern University. */

/**
 * \file simulate_model.cc
 * \brief Implements the molstat::SimulateModel class for simulating
 *    histograms.
 *
 * \author Matthew G.\ Reuter
 * \date October 2014
 */

#include "simulate_model.h"
#include "simulator_exceptions.h"

namespace molstat {

SimulateModelType SimulateModel::getModelType() const
{
	return std::type_index{ typeid(SimulateModel) };
}

std::size_t SimulateModel::get_num_parameters() const
{
	return get_names().size();
}

ObservableFunction SimulateModel::getObservableFunction(
	const ObservableIndex &obs) const
{
	ObservableFunction obsfunc;

	try
	{
		// get the observable's factory and then bind it to this
		obsfunc = (compatible_observables.at(obs))(shared_from_this());
	}
	catch(const std::out_of_range &e)
	{
		throw IncompatibleObservable();
	}

	return obsfunc;
}

std::valarray<double> SimulateModel::generateParameters(Engine &engine) const
{
	const std::size_t length = get_num_parameters();
	std::valarray<double> ret(length);

	for(std::size_t j = 0; j < length; ++j)
	{
		ret[j] = dists[j]->sample(engine);
	}

	return ret;
}

} // namespace molstat
