/* This file is a part of MolStat, which is distributed under the Creative
   Commons Attribution-NonCommercial 4.0 International Public License.
   MolStat (c) 2014, Northwestern University. */
/**
 * \file transport_observables.h
 * \brief Interfaces for the various observables related to electron transport.
 *
 * Each observable should have a function that returns either double or
 * std::array<double, N> after taking in a GSL rng handle and a vector<double>.
 *
 * \author Matthew G.\ Reuter
 * \date September 2014
 */

#ifndef __observables_h__
#define __observables_h__

#include <array>
#include <general/simulator_tools/simulate_model_interface.h>

namespace molstat {
namespace transport {

/**
 * \brief Observable class for the applied bias.
 *
 * \tparam MPs The number of model parameters used to calculates the applied
 *    bias.
 */
template<std::size_t MPs>
class AppliedBias : public Observable<MPs> {
public:
	AppliedBias() = default;
	virtual ~AppliedBias() = default;

	/**
	 * \brief Returns the applied bias for a set of model parameters.
	 *
	 * \param[in] params A set of model parameters.
	 * \return The applied bias for the model parameters.
	 */
	virtual double AppBias(const std::array<double, MPs> &params) const = 0;

	/**
	 * \brief Route the operator() to AppBias.
	 *
	 * \param[in] params A set of model parameters.
	 * \return The applied bias for the model parameters.
	 */
	virtual double operator()(const std::array<double, MPs> &params) const
		override final {

		return AppBias(params);
	}
};

/**
 * \brief Observable class for the static conductance.
 *
 * \tparam MPs The number of model parameters used to calculate the static
 *    conductance.
 */
template<std::size_t MPs>
class StaticConductance : public Observable<MPs> {
public:
	StaticConductance() = default;
	virtual ~StaticConductance() = default;

	/**
	 * \brief Returns the static conductance for a set of model parameters.
	 *
	 * \param[in] params A set of model parameters.
	 * \return The static conductance for the model parameters.
	 */
	virtual double StaticG(const std::array<double, MPs> &params) const = 0;

	/**
	 * \brief Route the operator() to StaticG.
	 *
	 * \param[in] params A set of model parameters.
	 * \return The static conductance for the model parameters.
	 */
	virtual double operator()(const std::array<double, MPs> &params) const
		override final {

		return StaticG(params);
	}
};

/**
 * \brief Observable class for the differential conductance.
 *
 * \tparam MPs The number of model parameters used to calculated the
 *    differnetial conductance.
 */
template<std::size_t MPs>
class DifferentialConductance : public Observable<MPs> {
public:
	DifferentialConductance() = default;
	virtual ~DifferentialConductance() = default;

	/**
	 * \brief Returns the static conductance for a set of model parameters.
	 *
	 * \param[in] params A set of model parameters.
	 * \return The differential conductance for the model parameters.
	 */
	virtual double DiffG(const std::array<double, MPs> &params) const = 0;

	/**
	 * \brief Route the operator() to DiffG.
	 *
	 * \param[in] params A set of model parameters.
	 * \return The differential conductance for the model parameters.
	 */
	virtual double operator()(const std::array<double, MPs> &params) const
		override final {

		return DiffG(params);
	}
};

} // namespace molstat::transport
} // namespace molstat

#endif
