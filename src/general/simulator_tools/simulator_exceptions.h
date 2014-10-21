/* This file is a part of MolStat, which is distributed under the Creative
   Commons Attribution-NonCommercial 4.0 International Public License.
   MolStat (c) 2014, Northwestern University. */
/**
 * \file simulator_exceptions.h
 * \brief Defines exceptions used by the simulator program.
 *
 * \author Matthew G.\ Reuter
 * \date October 2014
 */

#ifndef __simulator_exceptions_h__
#define __simulator_exceptions_h__

#include <stdexcept>
#include <string>

namespace molstat {

/**
 * \brief Exception thrown when a model is missing one of its required
 *    molstat::RandomDistributions.
 */
class MissingDistribution : public std::logic_error {
public:
	MissingDistribution() = delete;
	virtual ~MissingDistribution() = default;

	/**
	 * \brief Constructor that creates the error message, given the name of
	 *    the distribution.
	 *
	 * \param[in] name The name of the missing distribution.
	 */
	MissingDistribution(const std::string &name)
		: std::logic_error("Missing the distribution \"" + name + "\".") {
	}
};

/**
 * \brief Exception thrown when the desired observable is incompatible with
 *    the model used to simulate data.
 */
class IncompatibleObservable : public std::logic_error {
public:
	IncompatibleObservable() :
		std::logic_error("Incompatible observable.") {
	}
};

/**
 * \brief Exception thrown when no observables have been requested.
 */
class NoObservables : public std::logic_error {
public:
	NoObservables() :
		std::logic_error("No observables specified.") {
	}
};

/**
 * \brief Exception thrown when a molstat::CompositeObservable is requested
 *    by a model that is not a molstat::CompositeSimulateModel.
 */
class NotCompositeSimulateModel : public std::logic_error {
public:
	NotCompositeSimulateModel() :
		std::logic_error("Not a composite model.") {
	}
};

/**
 * \brief Exception thrown when a molstat::SimulateModel is not the correct type
 *    to be used with the specified molstat::CompositeSimulateModel.
 */
class IncompatibleSubmodel : public std::logic_error {
public:
	IncompatibleSubmodel() :
		std::logic_error("Incompatible submodel.") {
	}
};

} // namespace molstat

#endif