/* This file is a part of MolStat, which is distributed under the Creative
   Commons Attribution-NonCommercial 4.0 International Public License.

   (c) 2014 Northwestern University. */

/**
 * \file sym_two_site_channel.h
 * \brief Tight-binding channel with a two-site chain that couples symmetrically
 *    to both electrodes. The chain does not drop voltage.
 *
 * \author Matthew G.\ Reuter
 * \date October 2014
 */

#ifndef __sym_two_site_channel_h__
#define __sym_two_site_channel_h__

#include "observables.h"
#include "junction.h"

namespace molstat {
namespace transport {

/**
 * \brief Simulator submodel for transport through a two-site site that couples
 *    symmetrically to both electrodes.
 *
 * Inherited model parameters (from molstat::transport::TransportJunction) are
 * - `ef` (\f$E_\mathrm{F}\f$), the Fermi energy,
 * - `v` (\f$V\f$), the applied bias.
 *
 * Model parameters are
 * - `epsilon` (\f$\varepsilon\f$), the site-energy,
 * - `gamma` (\f$\Gamma\f$), the site/lead coupling,
 * - `beta` (\f$\beta\f$), the inter-site coupling.
 *
 * Starting from \f[ \hat{H} = \left[ \begin{array}{cc} \varepsilon & \beta \\ \beta & \varepsilon \end{array} \right], \;\;\;\;
 * \hat{\Sigma}_\mathrm{L} = \left[ \begin{array}{cc} -i\Gamma/2 & 0 \\ 0 & 0 \end{array} \right], \;\;\;\;
 * \hat{\Sigma}_\mathrm{R} = \left[ \begin{array}{cc} 0 & 0 \\ 0 & -i\Gamma/2 \end{array} \right], \f]
 * the transmission function is
 * \f[ T(E) = \frac{16 \Gamma^2 \beta^2}{\left[ 4(E-\varepsilon)^2-4\beta^2-\Gamma^2\right]^2 + 16 \Gamma^2(E-\varepsilon)^2}. \f]
 * - Electric current:
 *   \f{eqnarray*}{ G_\mathrm{s}(V) & = & \frac{2e}{h} \frac{2\beta\Gamma}{4\beta^2+\Gamma^2} \mathrm{Re} \left[ (\Gamma + i2\beta) \mathrm{arctanh}\left( \frac{2(E_\mathrm{F} - \varepsilon + eV/2)}{2\beta + i\Gamma} \right) \right] \\
 *   && -\frac{2e}{h} \frac{2\beta\Gamma}{4\beta^2+\Gamma^2} \mathrm{Re} \left[ (\Gamma + i2\beta) \mathrm{arctanh}\left( \frac{2(E_\mathrm{F} - \varepsilon - eV/2)}{2\beta + i\Gamma} \right) \right]. \f}
 * - Differential conductance:
 *   \f[ G_\mathrm{d}(V) = \frac{2e^2}{h} \frac{1}{2} \left[ T(E_\mathrm{F} + eV/2) + T(E_\mathrm{F} - eV/2) \right]. \f]
 * - Static conductance:
 *   \f{eqnarray*}{ G_\mathrm{s}(V) & = & \frac{2e^2}{h} \frac{2\beta\Gamma}{eV(4\beta^2+\Gamma^2)} \mathrm{Re} \left[ (\Gamma + i2\beta) \mathrm{arctanh}\left( \frac{2(E_\mathrm{F} - \varepsilon + eV/2)}{2\beta + i\Gamma} \right) \right] \\
 *   && -\frac{2e^2}{h} \frac{2\beta\Gamma}{eV(4\beta^2+\Gamma^2)} \mathrm{Re} \left[ (\Gamma + i2\beta) \mathrm{arctanh}\left( \frac{2(E_\mathrm{F} - \varepsilon - eV/2)}{2\beta + i\Gamma} \right) \right]. \f}
 */
class SymTwoSiteChannel : public Channel,
	public ElectricCurrent,
	public ZeroBiasConductance,
	public DifferentialConductance,
	public StaticConductance,
	public SeebeckCoefficient
{
private:
	/**
	 * \brief Calculates the antiderivative needed for the electric current
	 *    (fixed values of the model parameters).
	 *
	 * \param[in] z The limit of integration.
	 * \param[in] eps The channel energy, epsilon.
	 * \param[in] gamma The channel-lead coupling, gamma.
	 * \param[in] beta The site-site coupling, beta.
	 * \return The antiderivative needed for the static conductance.
	 */
	static double current_integral(const double z, const double eps,
		const double gamma, const double beta);

public:
	/// Container index for the Fermi energy.
	static const std::size_t Index_EF;

	/// Container index for the applied bias.
	static const std::size_t Index_V;

	/// Container index for the site energy.
	static const std::size_t Index_epsilon;

	/// Container index for the site-lead coupling.
	static const std::size_t Index_gamma;

	/// Container index for the inter-site coupling.
	static const std::size_t Index_beta;

protected:
	virtual std::vector<std::string> get_names() const override;

public:
	virtual ~SymTwoSiteChannel() = default;

	/**
	 * \brief Calculates the transmission for a set of model parameters.
	 *
	 * \param[in] e The energy of the incident electron.
	 * \param[in] v The applied bias.
	 * \param[in] eps The level energy.
	 * \param[in] gamma The level-electrode coupling.
	 * \param[in] beta The site-site coupling.
	 * \return The transmission for this set of parameters.
	 */
	static double transmission(const double e, const double v, const double eps,
		const double gamma, const double beta);

	virtual double ECurrent(const std::valarray<double> &params) const override;
	virtual double ZeroBiasG(const std::valarray<double> &params) const
		override;
	virtual double DiffG(const std::valarray<double> &params) const override;
	virtual double StaticG(const std::valarray<double> &params) const override;
	virtual double SeebeckS(const std::valarray<double> &params) const override;
};

} // namespace molstat::transport
} // namespace molstat

#endif
