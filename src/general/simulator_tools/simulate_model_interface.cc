/* This file is a part of MolStat, which is distributed under the Creative
   Commons Attribution-NonCommercial 4.0 International Public License.
   MolStat (c) 2014, Northwestern University. */
/**
 * \file simulate_model_interface.cc
 * \brief Implements aspects of abstract classes for simulating histograms.
 *
 * \author Matthew G.\ Reuter
 * \date October 2014
 */

// no includes because this file is included by simulate_model_interface.h

namespace molstat {

template<std::size_t MPs>
SimulateModel<MPs>::SimulateModel(
	const std::map<std::string,
	               std::shared_ptr<RandomDistribution>> &avail,
	const std::array<std::string, MPs> &names)
	: dists{ { nullptr } } {

	// go through each of the requested parameters
	for(std::size_t j = 0; j < MPs; ++j) {
		try {
			dists[j] = avail.at(names[j]);
		}
		catch(const std::out_of_range &e) {
			throw std::runtime_error(names[j].c_str());
		}
	}
}

template<std::size_t MPs>
std::array<double, MPs>
	SimulateModel<MPs>::generateParameters(gsl_rng_ptr &r) const {

	std::array<double, MPs> ret;

	for(std::size_t j = 0; j < MPs; ++j)
		ret[j] = dists[j]->sample(r);

	return ret;
}

template<std::size_t OBS, std::size_t MPs>
std::array<double, OBS>
	SimulateObservables<OBS, MPs>::simulate(gsl_rng_ptr &r) const {

	std::array<double, OBS> ret;

	// generate the model parameters
	const std::array<double, MPs> params = model->generateParameters(r);

	// calculate each observable
	for(std::size_t j = 0; j < OBS; ++j)
		ret[j] = observables[j](params);

	return ret;
}

} // namespace molstat
