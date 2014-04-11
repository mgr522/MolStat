/**
 * \file symmetric_voltage_independent.h
 * \brief The symmetric-coupling, voltage-independent tight-binding model for
 *        calculating conductances.
 *
 * \author Matthew G.\ Reuter
 * \date April 2014
 */

#ifndef __symmetric_voltage_independent_h__
#define __symmetric_voltage_independent_h__

#include <memory>
#include "rng.h"
#include "model_interface.h"

using std::shared_ptr;

/**
 * \brief Class encapsulating the voltage-independent model (symmetric
 *    coupling).
 *
 * Note: the tight-binding model, that is, the transmission function, is
 * indepedent of the voltage; the conductances, however, are not.
 */
class SymmetricVoltageIndependentModel : public ConductanceModel {
protected:
	/**
	 * \brief Random distribution for epsilon, the channel energy.
	 */
	shared_ptr<const RandomDistribution> dist_eps;

	/**
	 * \brief Random distribution for gamma, the channel-lead coupling.
	 */
	shared_ptr<const RandomDistribution> dist_gamma;

public:
	/**
	 * \brief Function that reads in probability distributions from the input
	 *    stream and constructs a VoltageIndependentModel using the desired
	 *    parameters.
	 *
	 * \exception std::runtime_error if the input is invalid or incomplete.
	 *
	 * \param[in,out] f Input stream.
	 * \return Shared pointer to the new ConductanceModel.
	 */
	static shared_ptr<ConductanceModel> create_model(FILE *f);

	SymmetricVoltageIndependentModel() = delete;

	/**
	 * \brief Constructor specifying the necessary random distributions.
	 *
	 * \param[in] eta The distribution for eta.
	 * \param[in] eps The distribution for epsilon.
	 * \param[in] gamma The distribution for gamma.
	 */
	SymmetricVoltageIndependentModel(shared_ptr<const RandomDistribution> eta,
		shared_ptr<const RandomDistribution> eps,
		shared_ptr<const RandomDistribution> gamma);

	/**
	 * \brief Destructor.
	 */
	virtual ~SymmetricVoltageIndependentModel() = default;

	/**
	 * \brief Gets the static conductance for a random set of model parameters.
	 *
	 * \param[in] r The handle for GSL random number generation.
	 * \param[in] EF The Fermi energy.
	 * \param[in] V The voltage.
	 * \return The static conductance.
	 */
	virtual double static_conductance(shared_ptr<gsl_rng> r,
		const double EF, const double V) const;

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
		const double EF, const double V) const;

	/**
	 * \brief Calculates the transmission for fixed values of epsilon and gamma.
	 *
	 * \param[in] E The incident energy of the electron.
	 * \param[in] eps The channel energy, epsilon.
	 * \param[in] gamma The channel-lead coupling, gamma.
	 * \return The transmission.
	 */
	static double transmission(const double E, const double eps,
		const double gamma);

	/**
	 * \brief Calculates the static conductance for fixed values of the model
	 *    parameters.
	 *
	 * \param[in] EF The Fermi energy.
	 * \param[in] V The voltage.
	 * \param[in] eta The relative voltage drops at the leads.
	 * \param[in] eps The channel energy, epsilon.
	 * \param[in] gamma The channel-lead coupling, gamma.
	 * \return The static conductance.
	 */
	static double static_conductance(const double EF, const double V,
		const double eta, const double eps, const double gamma);

	/**
	 * \brief Calculates the differential conductance for fixed values of the
	 *    model parameters.
	 *
	 * \param[in] EF The Fermi energy.
	 * \param[in] V The voltage.
	 * \param[in] eta The relative voltage drops at the leads.
	 * \param[in] eps The channel energy, epsilon.
	 * \param[in] gamma The channel-lead coupling, gamma.
	 * \return The differential conductance.
	 */
	static double diff_conductance(const double EF, const double V,
		const double eta, const double eps, const double gamma);
};

#endif
