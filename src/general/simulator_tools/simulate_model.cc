/* This file is a part of MolStat, which is distributed under the Creative
   Commons Attribution-NonCommercial 4.0 International Public License.
   MolStat (c) 2014, Northwestern University. */
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
#include <general/string_tools.h>

namespace molstat {

ObservableFunction SimulateModel::getObservableFunction(
		const std::type_index &obs) const {

	ObservableFunction obsfunc;

	try {
		obsfunc = (compatible_observables.at(obs))(shared_from_this());
	}
	catch(const std::out_of_range &e) {
		throw IncompatibleObservable();
	}

	return obsfunc;
}

std::valarray<double> SimulateModel::generateParameters(gsl_rng_ptr &r) const {
	// make sure all distributions have been specified
	if(!distributions_specified) {
		const std::vector<std::string> names = get_names();
		const std::size_t length = names.size();

		for(std::size_t j = 0; j < length; ++j) {
			try {
				if(dists.at(j) == nullptr)
					throw std::out_of_range(""); // enter the exception handler
			}
			catch(const std::out_of_range &e) {
				throw MissingDistribution(names[j]);
			}
		}
	}

	const std::size_t length = get_num_parameters();
	std::valarray<double> ret(length);

	for(std::size_t j = 0; j < length; ++j) {
		ret[j] = dists[j]->sample(r);
	}

	return ret;
}

void SimulateModel::setDistribution(const std::string &name,
	std::shared_ptr<const RandomDistribution> dist) {

	// find what position this distribution is at
	const std::vector<std::string> names { get_names() };
	const std::size_t length = names.size();
	std::size_t pos;

	for(pos = 0; pos < length; ++pos) {
		if(to_lower(name) == to_lower(names[pos]))
			break;
	}

	if(pos == length) // this distribution is not needed
		return;

	// make sure dists is long enough
	if(dists.size() < length)
		dists.resize(length);

	dists[pos] = dist;

	// update the distributions_specified flag
	distributions_specified = true;
	for(pos = 0; pos < length; ++pos)
		distributions_specified =
			distributions_specified && (dists[pos] != nullptr);
}

} // namespace molstat