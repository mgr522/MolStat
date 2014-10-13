/* This file is a part of MolStat, which is distributed under the Creative
   Commons Attribution-NonCommercial 4.0 International Public License.
   MolStat (c) 2014, Northwestern University. */
/**
 * \file sym_two_site_simulate_model.h
 * \brief Tight-binding model of a two-site chain that couples symmetrically to
 *    both electrodes. The chain does not drop voltage.
 *
 * \author Matthew G.\ Reuter
 * \date October 2014
 */

#ifndef __sym_two_site_simulate_model_h__
#define __sym_two_site_simulate_model_h__

#include <array>
#include <memory>
#include <string>
#include <map>
#include <general/simulator_tools/simulate_model_interface.h>
#include "transport_observables.h"

namespace molstat {
namespace transport {

/**
 * \brief Simulator model for transport through a two-site site that couples
 *    symmetrically to both electrodes.
 *
 * Model parameters are
 * - `ef` (\f$E_\mathrm{F}\f$), the Fermi energy,
 * - `v` (\f$V\f$), the applied bias,
 * - `epsilon` (\f$\varepsilon\f$), the site-energy,
 * - `gamma` (\f$\Gamma\f$), the site/lead coupling,
 * - `beta` (\f$\beta\f$), the inter-site coupling.
 *
 * Starting from \f[ \hat{H} = \left[ \begin{array}{cc} \varepsilon & \beta \\ \beta & \varepsilon \end{array} \right], \;\;\;\;
 * \hat{\Sigma}_\mathrm{L} = \left[ \begin{array}{cc} -i\Gamma/2 & 0 \\ 0 & 0 \end{array} \right], \;\;\;\;
 * \hat{\Sigma}_\mathrm{R} = \left[ \begin{array}{cc} 0 & 0 \\ 0 & -i\Gamma/2 \end{array} \right], \f]
 * the transmission function is
 * \f[ T(E) = \frac{16 \Gamma^2 \beta^2}{\left[ 4(E-\varepsilon)^2-4\beta^2-\Gamma^2\right]^2 + 16 \Gamma^2(E-\varepsilon)^2}. \f]
 * - Differential conductance:
 *   \f[ G_\mathrm{d}(V) = \frac{2e^2}{h} \frac{1}{2} \left[ T(E_\mathrm{F} + eV/2) + T(E_\mathrm{F} - eV/2) \right]. \f]
 * - Static conductance:
 *   \f{eqnarray*}{ G_\mathrm{s}(V) & = & \frac{2e^2}{h} \frac{2\beta\Gamma}{eV(4\beta^2+\Gamma^2)} \mathrm{Re} \left[ (\Gamma + i2\beta) \mathrm{arctanh}\left( \frac{2(E_\mathrm{F} - \varepsilon + eV/2)}{2\beta + i\Gamma} \right) \right] \\
 *   && -\frac{2e^2}{h} \frac{2\beta\Gamma}{eV(4\beta^2+\Gamma^2)} \mathrm{Re} \left[ (\Gamma + i2\beta) \mathrm{arctanh}\left( \frac{2(E_\mathrm{F} - \varepsilon - eV/2)}{2\beta + i\Gamma} \right) \right]. \f}
 */
class SymTwoSiteSimulateModel : public SimulateModel<5>,
	public AppliedBias<5>,
	public DifferentialConductance<5>,
	public StaticConductance<5> {

private:
	/**
	 * \internal
	 * \brief Calculates the antiderivative needed for the static
	 *    conductance (fixed values of the model parameters).
	 *
	 * \param[in] z The limit of integration.
	 * \param[in] eps The channel energy, epsilon.
	 * \param[in] gamma The channel-lead coupling, gamma.
	 * \param[in] beta The site-site coupling, beta.
	 * \return The antiderivative needed for the static conductance.
	 * \endinternal
	 */
	static double static_c_integral(const double z, const double eps,
		const double gamma, const double beta);

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
	 * \brief Container index for the site-lead coupling.
	 */
	static const std::size_t Index_gamma;

	/**
	 * \brief Container index for the inter-site coupling.
	 */
	static const std::size_t Index_beta;

	SymTwoSiteSimulateModel() = delete;
	virtual ~SymTwoSiteSimulateModel() = default;

	/**
	 * \brief Constructor specifying the required parameters.
	 *
	 * \param[in] avail The available distributions, keyed by name.
	 */
	SymTwoSiteSimulateModel(
		const std::map<std::string, std::shared_ptr<RandomDistribution>> &avail);

	/**
	 * \brief Returns the applied bias for a set of model parameters.
	 *
	 * \param[in] params A set of model parameters.
	 * \return The applied bias.
	 */
	virtual double AppBias(const std::array<double, 5> &params) const override;

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

	/**
	 * \brief Returns the differential conductance for a set of model
	 *    parameters.
	 * 
	 * \param[in] params A set of model parameters.
	 * \return The differential conductance.
	 */
	virtual double DiffG(const std::array<double, 5> &params) const override;

	/**
	 * \brief Returns the static conductance for a set of model parameters.
	 * 
	 * \param[in] params A set of model parameters.
	 * \return The static conductance.
	 */
	virtual double StaticG(const std::array<double, 5> &params) const override;
};

} // namespace molstat::transport
} // namespace molstat

#endif
