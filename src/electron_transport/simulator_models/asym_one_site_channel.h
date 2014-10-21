/* This file is a part of MolStat, which is distributed under the Creative
   Commons Attribution-NonCommercial 4.0 International Public License.
   MolStat (c) 2014, Northwestern University. */
/**
 * \file asym_one_site_channel.h
 * \brief Tight-binding channel with one site that couples asymmetrically to
 *    both electrodes.
 *
 * \author Matthew G.\ Reuter
 * \date October 2014
 */

#ifndef __asym_one_site_channel_h__
#define __asym_one_site_channel_h__

#include "observables.h"
#include "junction.h"

namespace molstat {
namespace transport {

/**
 * \brief Simulator submodel for transport through a single site that couples
 *    asymmetrically to both electrodes.
 *
 * Inherited model parameters (from molstat::transport::TransportJunction) are
 * - `ef` (\f$E_\mathrm{F}\f$), the Fermi energy,
 * - `v` (\f$V\f$), the applied bias.
 *
 * Model parameters are
 * - `epsilon` (\f$\varepsilon\f$), the site-energy,
 * - `gammal` (\f$\Gamma_\mathrm{L}\f$), the site/lead coupling for one electrode,
 * - `gammar` (\f$\Gamma_\mathrm{R}\f$), the site/lead coupling for the other electrode,
 * - `a` (\f$a\f$), the scaling factor for the voltage-dependence of \f$\varepsilon\f$.
 *
 * Starting from \f$ \hat{H} = \left[ \varepsilon-aeV \right] \f$ and
 * \f$ \hat{\Sigma}_\mathrm{L/R} = \left[ -i\Gamma_\mathrm{L,R}/2 \right] \f$,
 * the transmission function is
 * \f[ T(E) = \frac{4\Gamma_\mathrm{L} \Gamma_\mathrm{R}}{4(E-\varepsilon-aeV)^2 + (\Gamma_\mathrm{L} + \Gamma_\mathrm{R})^2}. \f]
 * - Differential conductance:
 *   \f[ G_\mathrm{d}(V) = \frac{2e^2}{h} \left[ (1/2-a) T(E_\mathrm{F} + eV/2) + (1/2+a) T(E_\mathrm{F} - eV/2) \right]. \f]
 * - Static conductance:
 *   \f[ G_\mathrm{s}(V) = \frac{2e^2}{h} \frac{2\Gamma_\mathrm{L} \Gamma_\mathrm{R}}{eV(\Gamma_\mathrm{L} +\Gamma_\mathrm{R})} \left[ \arctan\left( \frac{2[E_\mathrm{F} - \varepsilon + (1/2-a) eV]}{\Gamma_\mathrm{L} + \Gamma_\mathrm{R}} \right) - \arctan\left( \frac{2[E_\mathrm{F} - \varepsilon - (1/2+a)eV]}{\Gamma_\mathrm{L} + \Gamma_\mathrm{R}} \right) \right]. \f]
 */
class AsymOneSiteChannel : public Channel,
	public DifferentialConductance,
	public StaticConductance {

public:
	/**
	 * \brief Container index for the Fermi energy.
	 */
	static const std::size_t Index_EF;

	/**
	 * \brief Container index for the applied bias.
	 */
	static const std::size_t Index_V;
	
	/**
	 * \brief Container index for the site energy.
	 */
	static const std::size_t Index_epsilon;
	
	/**
	 * \brief Container index for the left site-lead coupling.
	 */
	static const std::size_t Index_gammaL;
	
	/**
	 * \brief Container index for the right site-lead coupling.
	 */
	static const std::size_t Index_gammaR;
	
	/**
	 * \brief Container index for the bias drop scaling factor.
	 */
	static const std::size_t Index_a;

protected:
	/**
	 * \brief Gets the names of model parameters for this channel.
	 *
	 * \return The names of model parameters.
	 */
	virtual std::vector<std::string> get_names() const override;

public:
	virtual ~AsymOneSiteChannel() = default;

	/**
	 * \brief Calculates the transmission for a set of model parameters.
	 *
	 * \param[in] e The energy of the incident electron.
	 * \param[in] v The applied bias.
	 * \param[in] eps The level energy.
	 * \param[in] gammal The left level-electrode coupling.
	 * \param[in] gammar The right level-electrode coupling.
	 * \param[in] a The voltage drop scaling factor.
	 * \return The transmission for this set of parameters.
	 */
	static double transmission(const double e, const double v, const double eps,
		const double gammal, const double gammar, const double a);

	/**
	 * \brief Returns the differential conductance for a set of model
	 *    parameters.
	 * 
	 * \param[in] params A set of model parameters.
	 * \return The differential conductance.
	 */
	virtual double DiffG(const std::valarray<double> &params) const override;

	/**
	 * \brief Returns the static conductance for a set of model parameters.
	 * 
	 * \param[in] params A set of model parameters.
	 * \return The static conductance.
	 */
	virtual double StaticG(const std::valarray<double> &params) const override;
};

} // namespace molstat::transport
} // namespace molstat

#endif
