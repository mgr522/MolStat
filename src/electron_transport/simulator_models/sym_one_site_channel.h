/* This file is a part of MolStat, which is distributed under the Creative
   Commons Attribution-NonCommercial 4.0 International Public License.

   (c) 2014 Northwestern University. */

/**
 * \file sym_one_site_channel.h
 * \brief Tight-binding channel with one site that couples symmetrically to
 *    both electrodes.
 *
 * \author Matthew G.\ Reuter
 * \date October 2014
 */

#ifndef __sym_one_site_channel_h__
#define __sym_one_site_channel_h__

#include "observables.h"
#include "junction.h"

namespace molstat {
namespace transport {

/**
 * \brief Simulator submodel (molstat::transport::Channel) for transport
 *    through a single site that couples symmetrically to both electrodes.
 *
 * Inherited model parameters (from molstat::transport::TransportJunction) are
 * - `ef` (\f$E_\mathrm{F}\f$), the Fermi energy,
 * - `v` (\f$V\f$), the applied bias.
 *
 * Submodel parameters are
 * - `epsilon` (\f$\varepsilon\f$), the site-energy,
 * - `gamma` (\f$\Gamma\f$), the site/lead coupling,
 * - `a` (\f$a\f$), the scaling factor for the voltage drop.
 *
 * Starting from \f$ \hat{H} = \left[ \varepsilon + aeV \right] \f$ and
 * \f$ \hat{\Sigma}_\mathrm{L/R} = \left[ -i\Gamma/2 \right] \f$, the
 * transmission function is
 * \f[ T(E) = \frac{\Gamma^2}{(E-\varepsilon-aeV)^2 + \Gamma^2}. \f]
 * - Electric current:
 *   \f[ I(V) = \frac{2e}{h} \Gamma \left[ \arctan\left( \frac{E_\mathrm{F} - \varepsilon + (1/2-a) eV}{\Gamma} \right) - \arctan\left( \frac{E_\mathrm{F} - \varepsilon - (1/2+a) eV}{\Gamma} \right) \right]. \f]
 * - Differential conductance (with intermediate quantities):
 *   \f[ \frac{\partial}{\partial V}T(E) = \frac{2ea\Gamma^2(E - \varepsilon-aeV)}{[(E-\varepsilon-aeV)^2+\Gamma^2]^2}; \f]
 *   \f[ \frac{2e}{h} \int\limits_{E_\mathrm{F}-eV/2}^{E_\mathrm{F}+eV/2} \mathrm{d}E \frac{\partial}{\partial V} T(E) = \frac{2e^2a}{h} T(E_\mathrm{F}-eV/2) - \frac{2e^2a}{h}T(E_\mathrm{F}+eV/2); \f]
 *   \f[ G_\mathrm{d}(V) = \frac{2e^2}{h} \left[ (1/2-a) T(E_\mathrm{F} + eV/2) + (1/2+a) T(E_\mathrm{F} - eV/2) \right]. \f]
 * - Static conductance:
 *   \f[ G_\mathrm{s}(V) = \frac{2e^2}{h} \frac{\Gamma}{eV} \left[ \arctan\left( \frac{E_\mathrm{F} - \varepsilon + (1/2-a) eV}{\Gamma} \right) - \arctan\left( \frac{E_\mathrm{F} - \varepsilon - (1/2+a) eV}{\Gamma} \right) \right]. \f]
 * - Seebeck coefficient (zero bias):
 *   \f[ S = \frac{-2 (E_\mathrm{F}-\varepsilon)}{(E_\mathrm{F}-\varepsilon)^2 + \Gamma^2}. \f]
 */
class SymOneSiteChannel : public Channel,
	public ElectricCurrent,
	public ZeroBiasConductance,
	public DifferentialConductance,
	public StaticConductance,
	public SeebeckCoefficient
{
public:
	/// Container index for the Fermi energy.
	static const std::size_t Index_EF;

	/// Container index for the applied bias.
	static const std::size_t Index_V;

	/// Container index for the site energy.
	static const std::size_t Index_epsilon;

	/// Container index for the site-lead coupling.
	static const std::size_t Index_gamma;

	/// Container index for the bias drop scaling factor.
	static const std::size_t Index_a;

protected:
	virtual std::vector<std::string> get_names() const override;

public:
	virtual ~SymOneSiteChannel() = default;

	/**
	 * \brief Calculates the transmission for a set of model parameters.
	 *
	 * \param[in] e The energy of the incident electron.
	 * \param[in] V The applied bias.
	 * \param[in] eps The level energy.
	 * \param[in] gamma The level-electrode coupling.
	 * \param[in] a The voltage drop scaling factor.
	 * \return The transmission for this set of parameters.
	 */
	static double transmission(const double e, const double V, const double eps,
		const double gamma, const double a);
	
	virtual double ECurrent(const std::valarray<double> &params) const override;
	virtual double ZeroBiasG(const std::valarray<double> &params) const override;
	virtual double DiffG(const std::valarray<double> &params) const override;
	virtual double StaticG(const std::valarray<double> &params) const override;
	virtual double SeebeckS(const std::valarray<double> &params) const override;
};

} // namespace molstat::transport
} // namespace molstat

#endif
