/* This file is a part of MolStat, which is distributed under the Creative
   Commons Attribution-NonCommercial 4.0 International Public License.
   MolStat (c) 2014, Northwestern University. */
/**
 * \file electron_transport/simulator_models/observables.h
 * \brief Interfaces for the various observables related to electron transport.
 *
 * \author Matthew G.\ Reuter
 * \date October 2014
 */

#ifndef __transport_observables_h__
#define __transport_observables_h__

#include <valarray>
#include <general/simulator_tools/observable.h>

namespace molstat {
namespace transport {

/// Observable class for the applied bias.
class AppliedBias : public Observable<AppliedBias>
{
public:
	AppliedBias()
		: Observable<AppliedBias>(&AppliedBias::AppBias)
	{}

	virtual ~AppliedBias() = default;

	/**
	 * \brief Returns the applied bias for a set of model parameters.
	 *
	 * \param[in] params A set of model parameters.
	 * \return The applied bias for the model parameters.
	 */
	virtual double AppBias(const std::valarray<double> &params) const = 0;
};

/// Observable class for the electric current.
class ElectricCurrent : public Observable<ElectricCurrent>
{
public:
	ElectricCurrent()
		: Observable<ElectricCurrent>(&ElectricCurrent::ECurrent)
	{}

	virtual ~ElectricCurrent() = default;

	/**
	 * \brief Returns the electric current for a set of model parameters.
	 *
	 * \param[in] params A set of model parameters.
	 * \return The electric current for the model parameters.
	 */
	virtual double ECurrent(const std::valarray<double> &params) const = 0;
};

/// Observable class for the static conductance.
class StaticConductance : public Observable<StaticConductance>
{
public:
	StaticConductance()
		: Observable<StaticConductance>(&StaticConductance::StaticG)
	{}

	virtual ~StaticConductance() = default;

	/**
	 * \brief Returns the static conductance for a set of model parameters.
	 *
	 * \param[in] params A set of model parameters.
	 * \return The static conductance for the model parameters.
	 */
	virtual double StaticG(const std::valarray<double> &params) const = 0;
};

/// Observable class for the zero-bias conductance.
class ZeroBiasConductance : public Observable<ZeroBiasConductance>
{
public:
	ZeroBiasConductance()
		: Observable<ZeroBiasConductance>(&ZeroBiasConductance::ZeroBiasG)
	{}

	virtual ~ZeroBiasConductance() = default;

	/**
	 * \brief Returns the zero-bias conductance for a set of model parameters.
	 *
	 * \param[in] params A set of model parameters.
	 * \return The zero-bias conductance for the model parameters.
	 */
	virtual double ZeroBiasG(const std::valarray<double> &params) const = 0;
};

/**
 * \brief Observable class for the differential conductance.
 *
 * This observable supercedes the molstat::transport::ZeroBiasConductance by
 * also providing a bias-dependent conductance. The differential conductance,
 * evaluated at an applied bias of 0, should be identical to the zero-bias
 * conductance.
 */
class DifferentialConductance : public Observable<DifferentialConductance>
{
public:
	DifferentialConductance()
		: Observable<DifferentialConductance>(&DifferentialConductance::DiffG)
	{}

	virtual ~DifferentialConductance() = default;

	/**
	 * \brief Returns the differential conductance for a set of model parameters.
	 *
	 * \param[in] params A set of model parameters.
	 * \return The differential conductance for the model parameters.
	 */
	virtual double DiffG(const std::valarray<double> &params) const = 0;
};

} // namespace molstat::transport
} // namespace molstat

#endif
