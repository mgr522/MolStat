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
protected:
	/**
	 * \brief Random distribution for eta, the relative voltage drop.
	 */
	shared_ptr<const RandomDistribution> dist_eta;

public:
	ConductanceModel() = delete;

	/**
	 * \brief Constructor specifying the distribution for eta.
	 *
	 * \param[in] eta The distribution for eta.
	 */
	ConductanceModel(shared_ptr<const RandomDistribution> eta);

	/**
	 * \brief Destructor.
	 */
	virtual ~ConductanceModel() = default;

	/**
	 * \brief Gets the static conductance for a random set of model parameters.
	 *
	 * \param[in] r The handle for GSL random number generation.
	 * \param[in] EF The Fermi energy.
	 * \param[in] V The voltage.
	 * \return The static conductance.
	 */
	virtual double static_conductance(shared_ptr<gsl_rng> r,
		const double EF, const double V) const = 0;

	/**
	 * \brief Gets the differential conductance for a random set of model
	 *        parameters.
	 *
	 * \param[in] r The handle for GSL random number generation.
	 * \param[in] EF The Fermi energy.
	 * \param[in] V The voltage.
	 * \return The differential conductance.
	 */
	virtual double diff_conductance(shared_ptr<gsl_rng> r,
		const double EF, const double V) const = 0;
};

#endif
