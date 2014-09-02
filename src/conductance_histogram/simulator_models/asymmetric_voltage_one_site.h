/* This file is a part of MolStat, which is distributed under the Creative
   Commons Attribution-NonCommercial 4.0 International Public License.
   MolStat (c) 2014, Northwestern University. */
/**
 * \file asymmetric_voltage_one_site.h
 * \brief The voltage-dependent asymmetric-coupling, tight-binding model for
 *        calculating conductances.
 *
 * \author Matthew G.\ Reuter
 * \date September 2014
 */

#ifndef __asymmetric_voltage_one_site_h__
#define __asymmetric_voltage_one_site_h__

#include <memory>
#include <general/random_distributions/rng.h>
#include "conductance_model.h"

using std::shared_ptr;

/**
 * \brief Class encapsulating the voltage-dependent, asymmetric-coupling,
 *    one-site tight-binding model.
 *
 * Model parameters are
 * - `epsilon` (\f$\varepsilon\f$), the site-energy,
 * - `gammaL` (\f$\Gamma_\mathrm{L}\f$), the site/lead coupling for one electrode.
 * - `gammaR` (\f$\Gamma_\mathrm{R}\f$), the site/lead coupling for the other electrode.
 * - `a` (\f$a\f$), the scaling factor for the voltage-dependence of \f$\varepsilon\f$.
 *
 * Starting from \f$ \hat{H} = \left[ \varepsilon-aeV \right] \f$ and
 * \f$ \hat{\Sigma}_\mathrm{L/R} = \left[ -i\Gamma_\mathrm{L,R}/2 \right] \f$,
 * the transmission function is
 * \f[ T(E) = \frac{4\Gamma_\mathrm{L} \Gamma_\mathrm{R}}{4(E-\varepsilon-aeV)^2 + (\Gamma_\mathrm{L} + \Gamma_\mathrm{R})^2}. \f]
 * - Differential conductance:
 *   \f[ g(V) = \frac{2e^2}{h} \left[ (\eta-a) T(E_\mathrm{F} + \eta eV) + (a+1-\eta) T(E_\mathrm{F} + (\eta-1)eV) \right]. \f]
 * - Static conductance:
 *   \f[ g(V) = \frac{2e^2}{h} \frac{2\Gamma_\mathrm{L} \Gamma_\mathrm{R}}{eV(\Gamma_\mathrm{L} +\Gamma_\mathrm{R})} \left[ \arctan\left( \frac{2[E_\mathrm{F} - \varepsilon + (\eta-a) eV]}{\Gamma_\mathrm{L} + \Gamma_\mathrm{R}} \right) - \arctan\left( \frac{2[E_\mathrm{F} - \varepsilon + (\eta-1-a)eV]}{\Gamma_\mathrm{L} + \Gamma_\mathrm{R}} \right) \right]. \f]
 */
class AsymmetricVoltageOneSiteModel : public ConductanceModel {
protected:
	/**
	 * \internal
	 * \brief Random distribution for \f$\varepsilon\f$, the channel energy.
	 * \endinternal
	 */
	shared_ptr<const RandomDistribution> dist_eps;

	/**
	 * \internal
	 * \brief Random distribution for \f$\Gamma_\mathrm{L}\f$, one channel-lead
	 *    coupling.
	 * \endinternal
	 */
	shared_ptr<const RandomDistribution> dist_gammaL;

	/**
	 * \internal
	 * \brief Random distribution for \f$\Gamma_\mathrm{R}\f$, one channel-lead
	 *    coupling.
	 * \endinternal
	 */
	shared_ptr<const RandomDistribution> dist_gammaR;

	/**
	 * \internal
	 * \brief Random distribution for \f$a\f$, the strength of the voltage
	 *    dependence.
	 * \endinternal
	 */
	shared_ptr<const RandomDistribution> dist_a;

public:
	AsymmetricVoltageOneSiteModel() = delete;

	/**
	 * \brief Constructor requiring the random distributions.
	 *
	 * \param[in] eps The distribution for \f$\varepsilon\f$.
	 * \param[in] gammaL The distribution for \f$\Gamma_\mathrm{L}\f$.
	 * \param[in] gammaR The distribution for \f$\Gamma_\mathrm{R}\f$.
	 * \param[in] a The distribution for \f$a\f$.
	 */
	AsymmetricVoltageOneSiteModel(
		const shared_ptr<const RandomDistribution> &eps,
		const shared_ptr<const RandomDistribution> &gammaL,
		const shared_ptr<const RandomDistribution> &gammaR,
		const shared_ptr<const RandomDistribution> &a);

	/**
	 * \internal
	 * \brief Destructor.
	 * \endinternal
	 */
	virtual ~AsymmetricVoltageOneSiteModel() = default;

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
	 * \param[in] V The voltage.
	 * \param[in] eps The channel energy, \f$\varepsilon\f$.
	 * \param[in] gammaL The left channel-lead coupling,
	 *    \f$\Gamma_\mathrm{L}\f$.
	 * \param[in] gammaR The right channel-lead coupling,
	 *    \f$\Gamma_\mathrm{R}\f$.
	 * \param[in] a The strength of the voltage dependence, \f$a\f$.
	 * \return The transmission.
	 */
	static double transmission(const double E, const double V, const double eps,
		const double gammaL, const double gammaR, const double a);

	/**
	 * \brief Calculates the static conductance for fixed values of the model
	 *    parameters.
	 *
	 * \param[in] EF The Fermi energy.
	 * \param[in] V The voltage.
	 * \param[in] eta The relative voltage drops at the leads.
	 * \param[in] eps The channel energy, \f$\varepsilon\f$.
	 * \param[in] gammaL The left channel-lead coupling,
	 *    \f$\Gamma_\mathrm{L}\f$.
	 * \param[in] gammaR The right channel-lead coupling,
	 *    \f$\Gamma_\mathrm{R}\f$.
	 * \param[in] a The strength of the voltage dependence, \f$a\f$.
	 * \return The static conductance.
	 */
	static double static_conductance(const double EF, const double V,
		const double eta, const double eps, const double gammaL,
		const double gammaR, const double a);

	/**
	 * \brief Calculates the differential conductance for fixed values of the
	 *    model parameters.
	 *
	 * \param[in] EF The Fermi energy.
	 * \param[in] V The voltage.
	 * \param[in] eta The relative voltage drops at the leads.
	 * \param[in] eps The channel energy, \f$\varepsilon\f$.
	 * \param[in] gammaL The left channel-lead coupling,
	 *    \f$\Gamma_\mathrm{L}\f$.
	 * \param[in] gammaR The right channel-lead coupling,
	 *    \f$\Gamma_\mathrm{R}\f$.
	 * \param[in] a The strength of the voltage dependence, \f$a\f$.
	 * \return The differential conductance.
	 */
	static double diff_conductance(const double EF, const double V,
		const double eta, const double eps, const double gammaL,
		const double gammaR, const double a);
};

#endif
