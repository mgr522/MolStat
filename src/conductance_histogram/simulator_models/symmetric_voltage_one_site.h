/* This file is a part of MolStat, which is distributed under the Creative
   Commons Attribution-NonCommercial 4.0 International Public License.
   MolStat (c) 2014, Northwestern University. */
/**
 * \file symmetric_voltage_one_site.h
 * \brief The symmetric-coupling, voltage-dependent one-site tight-binding
 *        model for calculating conductances.
 *
 * \author Matthew G.\ Reuter
 * \date July 2014
 */

#ifndef __symmetric_voltage_one_site_h__
#define __symmetric_voltage_one_site_h__

#include <memory>
#include <general/random_distributions/rng.h>
#include "conductance_model.h"

using std::shared_ptr;

/**
 * \brief Class encapsulating the voltage-dependent, one-site model (symmetric
 *    coupling).
 *
 * Model parameters are
 * - `epsilon` (\f$\varepsilon\f$), the site-energy,
 * - `gamma` (\f$\Gamma\f$), the site/lead coupling.
 * - `a` (\f$a\f$), the strength of the voltage dependence.
 *
 * Starting from \f$ \hat{H} = \left[ \varepsilon + aeV \right] \f$ and
 * \f$ \hat{\Sigma}_\mathrm{L/R} = \left[ -i\Gamma/2 \right] \f$, the
 * transmission function is
 * \f[ T(E) = \frac{\Gamma^2}{(E-\varepsilon-aeV)^2 + \Gamma^2}. \f]
 * - Differential conductance (with intermediate quantities):
 *   \f[ \frac{\partial}{\partial V}T(E) = \frac{2ea\Gamma^2(E - \varepsilon-aeV)}{[(E-\varepsilon-aeV)^2+\Gamma^2]^2}; \f]
 *   \f[ \frac{2e}{h} \int\limits_{E_\mathrm{F}+(\eta-1)eV}^{E_\mathrm{F}+\eta eV} \mathrm{d}E \frac{\partial}{\partial V} T(E) = \frac{2e^2a}{h} T(E_\mathrm{F}+(\eta-1)eV) - \frac{2e^2a}{h}T(E_\mathrm{F}+\eta eV); \f]
 *   \f[ g(V) = \frac{2e^2}{h} \left[ (\eta-a) T(E_\mathrm{F} + \eta eV) + (a+1-\eta) T(E_\mathrm{F} + (\eta-1)eV) \right]. \f]
 * - Static conductance:
 *   \f[ g(V) = \frac{2e^2}{h} \frac{\Gamma}{eV} \left[ \arctan\left( \frac{E_\mathrm{F} - \varepsilon + (\eta-a) eV}{\Gamma} \right) - \arctan\left( \frac{E_\mathrm{F} - \varepsilon + (\eta-1-a) eV}{\Gamma} \right) \right]. \f]
 */
class SymmetricVoltageOneSiteModel : public ConductanceModel {
protected:
	/**
	 * \internal
	 * \brief Random distribution for \f$\varepsilon\f$, the channel energy.
	 * \endinternal
	 */
	shared_ptr<const RandomDistribution> dist_eps;

	/**
	 * \internal
	 * \brief Random distribution for \f$\Gamma\f$, the channel-lead coupling.
	 * \endinternal
	 */
	shared_ptr<const RandomDistribution> dist_gamma;

	/**
	 * \internal
	 * \brief Random distribution for \f$a\f$, the strength of the voltage
	 *    dependence.
	 * \endinternal
	 */
	shared_ptr<const RandomDistribution> dist_a;

public:
	SymmetricVoltageOneSiteModel() = delete;

	/**
	 * \brief Constructor requiring the random distributions.
	 *
	 * \param[in] eps The distribution for \f$\varepsilon\f$.
	 * \param[in] gamma The distribution for \f$\Gamma\f$.
	 * \param[in] a The distribution for \f$a\f$.
	 */
	SymmetricVoltageOneSiteModel(
		const shared_ptr<const RandomDistribution> &eps,
		const shared_ptr<const RandomDistribution> &gamma,
		const shared_ptr<const RandomDistribution> &a);

	/**
	 * \internal
	 * \brief Destructor.
	 * \endinternal
	 */
	virtual ~SymmetricVoltageOneSiteModel() = default;

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
	 *    model parameters.
	 *
	 * \param[in] r The handle for GSL random number generation.
	 * \param[in] EF The Fermi energy.
	 * \return the zero-bias (differential) conductance.
	 */
	virtual double zero_bias_conductance(shared_ptr<gsl_rng> r,
		const double EF) const;

	/**
	 * \brief Calculates the transmission for fixed values of epsilon and gamma.
	 *
	 * \param[in] E The incident energy of the electron.
	 * \param[in] V The voltage.
	 * \param[in] eps The channel energy, \f$\varepsilon\f$.
	 * \param[in] gamma The channel-lead coupling, \f$\Gamma\f$.
	 * \param[in] a The strength of the voltage dependence, \f$a\f$.
	 * \return The transmission.
	 */
	static double transmission(const double E, const double V, const double eps,
		const double gamma, const double a);

	/**
	 * \brief Calculates the static conductance for fixed values of the model
	 *    parameters.
	 *
	 * \param[in] EF The Fermi energy.
	 * \param[in] V The voltage.
	 * \param[in] eta The relative voltage drops at the leads.
	 * \param[in] eps The channel energy, \f$\varepsilon\f$.
	 * \param[in] gamma The channel-lead coupling, \f$\Gamma\f$.
	 * \param[in] a The strength of the voltage-dependence, \f$a\f$.
	 * \return The static conductance.
	 */
	static double static_conductance(const double EF, const double V,
		const double eta, const double eps, const double gamma, const double a);

	/**
	 * \brief Calculates the differential conductance for fixed values of the
	 *    model parameters.
	 *
	 * \param[in] EF The Fermi energy.
	 * \param[in] V The voltage.
	 * \param[in] eta The relative voltage drops at the leads.
	 * \param[in] eps The channel energy, \f$\varepsilon\f$.
	 * \param[in] gamma The channel-lead coupling, \f$\Gamma\f$.
	 * \param[in] a The strength of the voltage-dependence, \f$a\f$.
	 * \return The differential conductance.
	 */
	static double diff_conductance(const double EF, const double V,
		const double eta, const double eps, const double gamma, const double a);
};

#endif
