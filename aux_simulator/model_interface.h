/**
 * \file model_interface.h
 * \brief Interface for the various tight-binding models and their transmission
 *        functions/conductances.
 *
 * \author Matthew G.\ Reuter
 * \date April 2014
 */

#ifndef __model_interface_h__
#define __model_interface_h__

#include <memory>
#include "rng.h"

using std::shared_ptr;

/**
 * \brief Abstract base class a conductance model.
 */
class ConductanceModel {
public:
	/**
	 * \brief Default constructor.
	 */
	ConductanceModel() = default;

	/**
	 * \brief Destructor.
	 */
	virtual ~ConductanceModel() = default;

	/**
	 * \brief Gets the static conductance for a random set of model parameters.
	 *
	 * \param[in] r The handle for GSL random number generation.
	 * \param[in] EF The Fermi energy.
	 * \param[in] eta The relative voltage drop.
	 * \param[in] V The voltage.
	 * \return The static conductance.
	 */
	virtual double static_conductance(shared_ptr<gsl_rng> r,
		const double EF, const double eta, const double V) const = 0;

	/**
	 * \brief Gets the differential conductance for a random set of model
	 *        parameters.
	 *
	 * \param[in] r The handle for GSL random number generation.
	 * \param[in] EF The Fermi energy.
	 * \param[in] eta The relative voltage drop.
	 * \param[in] V The voltage.
	 * \return The differential conductance.
	 */
	virtual double diff_conductance(shared_ptr<gsl_rng> r,
		const double EF, const double eta, const double V) const = 0;
};

#endif
