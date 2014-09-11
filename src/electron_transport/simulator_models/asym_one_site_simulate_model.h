/* This file is a part of MolStat, which is distributed under the Creative
   Commons Attribution-NonCommercial 4.0 International Public License.
   MolStat (c) 2014, Northwestern University. */
/**
 * \file asym_one_site_simulate_model.h
 * \brief Tight-binding model with one site that couples asymmetrically to
 *    both electrodes.
 *
 * \author Matthew G.\ Reuter
 * \date September 2014
 */

#ifndef __asym_one_site_simulate_model_h__
#define __asym_one_site_simulate_model_h__

#include <memory>
#include <string>
#include <vector>
#include <map>
#include <general/random_distributions/rng.h>
#include <general/simulator_tools/simulate_model_interface.h>
#include "transport_observables.h"

using std::shared_ptr;

/**
 * \brief Simulator model for transport through a single site that couples
 *    asymmetrically to both electrodes.
 *
 * Model parameters are
 * - `ef` (\f$E_\mathrm{F}\f$), the Fermi energy,
 * - `v` (\f$V\f$), the applied bias,
 * - `epsilon` (\f$\varepsilon\f$), the site-energy,
 * - `gammaL` (\f$\Gamma_\mathrm{L}\f$), the site/lead coupling for one electrode,
 * - `gammaR` (\f$\Gamma_\mathrm{R}\f$), the site/lead coupling for the other electrode,
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
class AsymOneSiteSimulateModel : public SimulateModel,
	public DifferentialConductance, public StaticConductance {

private:
	/**
	 * \brief Ordered list (vector) of the parameters needed for this model.
	 *
	 * If this order is changed, AsymOneSiteSimulateModel::unpack_parameters
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
	 * \param[out] gammal The left site-lead coupling.
	 * \param[out] gammar The right site-lead coupling.
	 * \param[out] a The voltage drop parameter.
	 */
	static void unpack_parameters(const std::vector<double> &vec, double &ef,
		double &v, double &epsilon, double &gammal, double &gammar, double &a);

public:
	AsymOneSiteSimulateModel() = delete;
	virtual ~AsymOneSiteSimulateModel() = default;

	/**
	 * \internal
	 * \brief Constructor specifying the required parameters.
	 *
	 * \param[in] avail The available distributions, keyed by name.
	 * \endinternal
	 */
	AsymOneSiteSimulateModel(
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
	 * \param[in] a The voltage drop scaling factor.
	 * \return The transmission for this set of parameters.
	 */
	static double transmission(const double e, const double v, const double eps,
		const double gammal, const double gammar, const double a);

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
