/* This file is a part of MolStat, which is distributed under the Creative
   Commons Attribution-NonCommercial 4.0 International Public License.
   MolStat (c) 2014, Northwestern University. */
/**
 * \file simulator.h
 * \brief Defines the molstat::Simulator class for simulating histograms.
 *
 * \todo Write a better description of how all these classes work together.
 *
 * \author Matthew G.\ Reuter
 * \date October 2014
 */

#ifndef __simulator_h__
#define __simulator_h__

#include <memory>
#include <valarray>
#include <vector>
#include <typeindex>
#include <general/random_distributions/rng.h>
#include "simulate_model.h"

namespace molstat {

/**
 * \brief Class for simulating data.
 *
 * This class contains the actual model used to simulate data and functions
 * for constructing the model, setting observables, etc.
 */
class Simulator {
private:
	/**
	 * \brief The actual model used to simulate data.
	 */
	std::shared_ptr<SimulateModel> model;

	/**
	 * \brief The functions that calculate observables, given a set of model
	 *    parameters.
	 */
	std::vector<ObservableFunction> obs_functions;

public:
	Simulator() = delete;
	~Simulator() = default;

	/**
	 * \brief Constructs molstat::Simulator using the specified model.
	 *
	 * \param[in] model_ The model to be used.
	 */
	Simulator(std::shared_ptr<SimulateModel> &model_);

	/**
	 * \brief Calculates the desired observables using the random number
	 *    generator.
	 *
	 * \throw molstat::NoObservables if no observables have been set.
	 * \throw molstat::MissingDistribution if any of the model's required
	 *    distributions is unspecified.
	 *
	 * \param[in] r Handle to the GSL random number generator.
	 * \return The simulated observables.
	 */
	std::valarray<double> simulate(gsl_rng_ptr &r) const;

	/**
	 * \brief Sets the `j`th observable for the simulator.
	 *
	 * First verify that the observable is compatible with the specified model.
	 * If compatible, get a function for calculating the observable and store
	 * it for use later.
	 *
	 * \throw molstat::IncompatibleObservable If our models is incompatible
	 *    with the desired observable.
	 * \throw out_of_range If `j` is out of range (not between 0 and
	 *    obs_functions.size()).
	 *
	 * \param[in] j The output index for this observable.
	 * \param[in] obs The identifier of the observable.
	 */
	void setObservable(std::size_t j, const ObservableIndex &obs);

	/**
	 * \brief Function for a \"zero\" observable.
	 *
	 * \todo Delete this function if not needed.
	 *
	 * \param[in] params A set of model parameters.
	 * \return The observable; in this case, 0.
	 */
	constexpr static double ZeroObs(const std::valarray<double> &params)
		noexcept {

		return 0.;
	}
};

} // namespace MolStat

#endif