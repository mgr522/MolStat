/* This file is a part of MolStat, which is distributed under the Creative
   Commons Attribution-NonCommercial 4.0 International Public License.
   MolStat (c) 2014, Northwestern University. */
/**
 * \file asym_two_site_simulate_model.h
 * \brief Tight-binding model of a two-site chain that couples asymmetrically
 *    to both electrodes. The chain does not drop voltage.
 *
 * \author Matthew G.\ Reuter
 * \date September 2014
 */

#ifndef __asym_two_site_simulate_model_h__
#define __asym_two_site_simulate_model_h__

#include <memory>
#include <string>
#include <vector>
#include <map>
#include <general/random_distributions/rng.h>
#include <general/simulator_tools/simulate_model_interface.h>
#include "transport_observables.h"

using std::shared_ptr;

/**
 * \brief Simulator model for transport through a two-site site that couples
 *    asymmetrically to both electrodes.
 *
 * Model parameters are
 * - `ef` (\f$E_\mathrm{F}\f$), the Fermi energy,
 * - `v` (\f$V\f$), the applied bias,
 * - `epsilon` (\f$\varepsilon\f$), the site-energy,
 * - `gammaL` (\f$\Gamma_\mathrm{L}\f$), the site/lead coupling for one
 *   electrode,
 * - `gammaR` (\f$\Gamma_\mathrm{R}\f$), the site/lead coupling for the other
 *   electrode,
 * - `beta` (\f$\beta\f$), the inter-site coupling.
 *
 * Starting from \f[ \hat{H} = \left[ \begin{array}{cc} \varepsilon & \beta \\ \beta & \varepsilon \end{array} \right], \;\;\;\;
 * \hat{\Sigma}_\mathrm{L} = \left[ \begin{array}{cc} -i\Gamma_\mathrm{L}/2 & 0 \\ 0 & 0 \end{array} \right], \;\;\;\;
 * \hat{\Sigma}_\mathrm{R} = \left[ \begin{array}{cc} 0 & 0 \\ 0 & -i\Gamma_\mathrm{R}/2 \end{array} \right], \f]
 * the transmission function is
 * \f[ T(E) = \frac{16 \Gamma_\mathrm{L} \Gamma_\mathrm{R} \beta^2}{\left[ 4(E-\varepsilon)^2-4\beta^2-\Gamma_\mathrm{L} \Gamma_\mathrm{R} \right]^2 + 4 (\Gamma_\mathrm{L} + \Gamma_\mathrm{R})^2(E-\varepsilon)^2}. \f]
 * - Differential conductance:
 *   \f[ G_\mathrm{d}(V) = \frac{2e^2}{h} \frac{1}{2} \left[ T(E_\mathrm{F} + eV/2) + T(E_\mathrm{F} - eV/2) \right]. \f]
 * - Indefinite integral for the static conductance:
 *   \f{eqnarray*}{ \int \mathrm{d}E T(E) & = & \frac{8\sqrt{2}\Gamma_\mathrm{L} \Gamma_\mathrm{R} \beta^2}{(\Gamma_\mathrm{L}+\Gamma_\mathrm{R})\sqrt{(\Gamma_\mathrm{L} - \Gamma_\mathrm{R})^2 - 16\beta^2}} \left\{ \frac{\mathrm{arctan}\left[ \frac{\sqrt{8}(E-\varepsilon)}{\sqrt{\Gamma_\mathrm{L}^2 + \Gamma_\mathrm{R}^2 - 8\beta^2 - (\Gamma_\mathrm{L} + \Gamma_\mathrm{R})\sqrt{(\Gamma_\mathrm{L} - \Gamma_\mathrm{R})^2 - 16\beta^2}}} \right]}{\sqrt{\Gamma_\mathrm{L}^2 + \Gamma_\mathrm{R}^2 -8\beta^2 - (\Gamma_\mathrm{L} + \Gamma_\mathrm{R}) \sqrt{(\Gamma_\mathrm{L} - \Gamma_\mathrm{R})^2-16\beta^2}}} \right. \\
 *   && \left. - \frac{\mathrm{arctan}\left[ \frac{\sqrt{8}(E-\varepsilon)}{\sqrt{\Gamma_\mathrm{L}^2 + \Gamma_\mathrm{R}^2 - 8\beta^2 + (\Gamma_\mathrm{L} + \Gamma_\mathrm{R})\sqrt{(\Gamma_\mathrm{L} - \Gamma_\mathrm{R})^2 - 16\beta^2}}} \right]}{\sqrt{\Gamma_\mathrm{L}^2 + \Gamma_\mathrm{R}^2 -8\beta^2 + (\Gamma_\mathrm{L} + \Gamma_\mathrm{R}) \sqrt{(\Gamma_\mathrm{L} - \Gamma_\mathrm{R})^2-16\beta^2}}} \right\}. \f}
 */
class AsymTwoSiteSimulateModel : public SimulateModel,
	public DifferentialConductance, public StaticConductance {

private:
	/**
	 * \brief Ordered list (vector) of the parameters needed for this model.
	 *
	 * If this order is changed, AsymTwoSiteSimulateModel::unpack_parameters
	 * also needs to be updated.
	 */
	static const std::vector<std::string> parameters;

	/**
	 * \brief Unpack a set of parameters from a vector to doubles.
	 *
	 * \param[in] vec The vector containing a set of parameters.
	 * \param[out] ef The Fermi energy.
	 * \param[out] v The voltage.
	 * \param[out] epsilon The site energy.
	 * \param[out] gammal The left lead-site coupling.
	 * \param[out] gammar The right lead-site coupling.
	 * \param[out] beta The site-site coupling.
	 */
	static void unpack_parameters(const std::vector<double> &vec, double &ef,
		double &v, double &epsilon, double &gammal, double &gammar,
		double &beta);

	/**
	 * \internal
	 * \brief Calculates the antiderivative needed for the static
	 *    conductance (fixed values of the model parameters).
	 *
	 * \param[in] z The limit of integration.
	 * \param[in] eps The channel energy, \f$\varepsilon\f$.
	 * \param[in] gammal The left channel-lead coupling,
	 *    \f$\Gamma_\mathrm{L}\f$.
	 * \param[in] gammar The right channel-lead coupling,
	 *    \f$\Gamma_\mathrm{R}\f$.
	 * \param[in] beta The site-site coupling, \f$\beta\f$.
	 * \return The antiderivative needed for the static conductance.
	 * \endinternal
	 */
	static double static_c_integral(const double z, const double eps,
		const double gammal, const double gammar, const double beta);

public:
	AsymTwoSiteSimulateModel() = delete;
	virtual ~AsymTwoSiteSimulateModel() = default;

	/**
	 * \internal
	 * \brief Constructor specifying the required parameters.
	 *
	 * \param[in] avail The available distributions, keyed by name.
	 * \endinternal
	 */
	AsymTwoSiteSimulateModel(
		const std::map<std::string, shared_ptr<RandomDistribution>> &avail);

	/**
	 * \brief Returns the differential conductance for a randomly-generated set
	 *    of model parameters.
	 * 
	 * \param[in] r The GSL random number generator handle.
	 * \return The differential conductance.
	 */
	virtual std::array<double, 2> DiffG(shared_ptr<gsl_rng> r) const
		override;

	/**
	 * \brief Returns the static conductance for a randomly-generated set of
	 *    model parameters.
	 * 
	 * \param[in] r The GSL random number generator handle.
	 * \return The static conductance.
	 */
	virtual std::array<double, 2> StaticG(shared_ptr<gsl_rng> r) const
		override;

	/**
	 * \brief Calculates the transmission for a set of model parameters.
	 *
	 * \param[in] e The energy of the incident electron.
	 * \param[in] v The applied bias.
	 * \param[in] eps The level energy.
	 * \param[in] gammal The left level-electrode coupling.
	 * \param[in] gammar The right level-electrode coupling.
	 * \param[in] beta The site-site coupling.
	 * \return The transmission for this set of parameters.
	 */
	static double transmission(const double e, const double v, const double eps,
		const double gammal, const double gammar, const double beta);

	/**
	 * \brief Calculates the static conductance for a set of model parameters.
	 *
	 * \param[in] vec The vector of model parameters.
	 * \return The static conductance for this set of parameters.
	 */
	static double static_conductance(const std::vector<double> &vec);

	/**
	 * \brief Calculates the differential conductance for a set of model
	 *    parameters.
	 *
	 * \param[in] vec The vector of model parameters.
	 * \return The differential conductance for this set of parameters.
	 */
	static double diff_conductance(const std::vector<double> &vec);
};

#endif
