/**
 * \file asymmetric_voltage_independent.h
 * \brief The asymmetric-coupling, voltage-independent tight-binding model for
 *        calculating conductances.
 *
 * \author Matthew G.\ Reuter
 * \date April 2014
 */

#ifndef __asymmetric_voltage_independent_h__
#define __asymmetric_voltage_independent_h__

#include <memory>
#include "rng.h"
#include "conductance_model.h"

using std::shared_ptr;

/**
 * \brief Class encapsulating the voltage-independent model (asymmetric
 *    coupling).
 *
 * Note: the tight-binding model, that is, the transmission function, is
 * indepedent of the voltage; the conductances, however, are not.
 */
class AsymmetricVoltageIndependentModel : public ConductanceModel {
protected:
	/**
	 * \brief Random distribution for epsilon, the channel energy.
	 */
	shared_ptr<const RandomDistribution> dist_eps;

	/**
	 * \brief Random distribution for gammaL, one channel-lead coupling.
	 */
	shared_ptr<const RandomDistribution> dist_gammaL;

	/**
	 * \brief Random distribution for gammaR, one channel-lead coupling.
	 */
	shared_ptr<const RandomDistribution> dist_gammaR;

public:
	AsymmetricVoltageIndependentModel() = delete;

	/**
	 * \brief Constructor requiring the random distributions.
	 *
	 * \param[in] eps The distribution for epsilon.
	 * \param[in] gammaL The distribution for gammaL.
	 * \param[in] gammaR The distribution for gammaR.
	 */
	AsymmetricVoltageIndependentModel(
		const shared_ptr<const RandomDistribution> &eps,
		const shared_ptr<const RandomDistribution> &gammaL,
		const shared_ptr<const RandomDistribution> &gammaR);

	/**
	 * \brief Destructor.
	 */
	virtual ~AsymmetricVoltageIndependentModel() = default;

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
		const double EF, const double eta, const double V) const;

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
		const double EF, const double eta, const double V) const;

	/**
	 * \brief Gets the zero-bias (differential) conductance for a random set of
	 *        model parameters.
	 *
	 * \param[in] r The handle for GSL random number generation.
	 * \param[in] EF The Fermi energy.
	 * \return The differential conductance.
	 */
	virtual double zero_bias_conductance(shared_ptr<gsl_rng> r,
		const double EF) const;

	/**
	 * \brief Calculates the transmission for fixed values of epsilon and the.
	 *     gammas.
	 *
	 * \param[in] E The incident energy of the electron.
	 * \param[in] eps The channel energy, epsilon.
	 * \param[in] gammaL The left channel-lead coupling, gammaL.
	 * \param[in] gammaR The right channel-lead coupling, gammaR.
	 * \return The transmission.
	 */
	static double transmission(const double E, const double eps,
		const double gammaL, const double gammaR);

	/**
	 * \brief Calculates the static conductance for fixed values of the model
	 *    parameters.
	 *
	 * \param[in] EF The Fermi energy.
	 * \param[in] V The voltage.
	 * \param[in] eta The relative voltage drops at the leads.
	 * \param[in] eps The channel energy, epsilon.
	 * \param[in] gammaL The left channel-lead coupling, gammaL.
	 * \param[in] gammaR The right channel-lead coupling, gammaR.
	 * \return The static conductance.
	 */
	static double static_conductance(const double EF, const double V,
		const double eta, const double eps, const double gammaL,
		const double gammaR);

	/**
	 * \brief Calculates the differential conductance for fixed values of the
	 *    model parameters.
	 *
	 * \param[in] EF The Fermi energy.
	 * \param[in] V The voltage.
	 * \param[in] eta The relative voltage drops at the leads.
	 * \param[in] eps The channel energy, epsilon.
	 * \param[in] gammaL The left channel-lead coupling, gammaL.
	 * \param[in] gammaR The right channel-lead coupling, gammaR.
	 * \return The differential conductance.
	 */
	static double diff_conductance(const double EF, const double V,
		const double eta, const double eps, const double gammaL,
		const double gammaR);
};

#endif
