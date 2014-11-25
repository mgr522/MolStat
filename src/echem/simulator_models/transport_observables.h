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

#include <memory>
#include <array>
#include <gsl/gsl_rng.h>
#include <general/simulator_tools/simulate_model_interface.h>

using std::shared_ptr;

/**
 * \brief Interface for the static conductance.
 */
class StaticConductance {
public:
	StaticConductance() = default;
	virtual ~StaticConductance() = default;

	/**
	 * \brief Returns the applied bias and static conductance for a randomly-
	 *    generated set of model parameters.
	 *
	 * \param[in] r The GSL random number generator handle.
	 * \return Array containing the applied bias (0) and the static conductance
	 *    (1).
	 */
	virtual std::array<double, 2> StaticG(shared_ptr<gsl_rng> r) const = 0;
};

/**
 * \brief Interface for the differential conductance.
 */
class DifferentialConductance {
public:
	DifferentialConductance() = default;
	virtual ~DifferentialConductance() = default;

	/**
	 * \brief Returns the applied bias and static conductance for a randomly-
	 *    generated set of model parameters.
	 *
	 * \param[in] r The GSL random number generator handle.
	 * \return Array containing the applied bias (0) and the differential
	 *    conductance (1).
	 */
	virtual std::array<double, 2> DiffG(shared_ptr<gsl_rng> r) const = 0;
};

/**
 * \brief Interface for the single molecule cyclic voltammetry current spike
 */
class SingMolCVPeak{
public:
    SingMolCVPeak() = default;
    virtual ~SingMolCVPeak() = default;

    /**
     * \brief Returns current peak potentials for a randomly-generated set of model parameters.
     *
     * \param[in] r the GSL random number generator handle.
     * \return Array containing the current peak potential during the forward scanning (0) and the backward scanning (1).
     */
    virtual std::array<double, 2> PeakPotentials(shared_ptr<gsl_rng> r) const = 0;
};
#endif
