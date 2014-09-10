/* This file is a part of MolStat, which is distributed under the Creative
   Commons Attribution-NonCommercial 4.0 International Public License.
   MolStat (c) 2014, Northwestern University. */
/**
 * \file sym_one_site_simulate_model.h
 * \brief Tight-binding model with one site that couples symmetrically to
 *    both electrodes.
 *
 * \author Matthew G.\ Reuter
 * \date September 2014
 */

#ifndef __sym_one_site_simulate_model_h__
#define __sym_one_site_simulate_model_h__

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
 *    symmetrically to both electrodes.
 *
 * Model parameters are
 * - `ef` (\f$E_\mathrm{F}\f$), the Fermi energy,
 * - `v` (\f$V\f$), the applied bias,
 * - `epsilon` (\f$\varepsilon\f$), the site-energy,
 * - `gamma` (\f$\Gamma\f$), the site/lead coupling,
 * - `a` (\f$a\f$), the scaling factor for the voltage drop.
 *
 * Starting from \f$ \hat{H} = \left[ \varepsilon + aeV \right] \f$ and
 * \f$ \hat{\Sigma}_\mathrm{L/R} = \left[ -i\Gamma/2 \right] \f$, the
 * transmission function is
 * \f[ T(E) = \frac{\Gamma^2}{(E-\varepsilon-aeV)^2 + \Gamma^2}. \f]
 * - Differential conductance (with intermediate quantities):
 *   \f[ \frac{\partial}{\partial V}T(E) = \frac{2ea\Gamma^2(E - \varepsilon-aeV)}{[(E-\varepsilon-aeV)^2+\Gamma^2]^2}; \f]
 *   \f[ \frac{2e}{h} \int\limits_{E_\mathrm{F}-eV/2}^{E_\mathrm{F}+eV/2} \mathrm{d}E \frac{\partial}{\partial V} T(E) = \frac{2e^2a}{h} T(E_\mathrm{F}-eV/2) - \frac{2e^2a}{h}T(E_\mathrm{F}+eV/2); \f]
 *   \f[ G_\mathrm{d}(V) = \frac{2e^2}{h} \left[ (1/2-a) T(E_\mathrm{F} + eV/2) + (1/2+a) T(E_\mathrm{F} - eV/2) \right]. \f]
 * - Static conductance:
 *   \f[ G_\mathrm{s}(V) = \frac{2e^2}{h} \frac{\Gamma}{eV} \left[ \arctan\left( \frac{E_\mathrm{F} - \varepsilon + (1/2-a) eV}{\Gamma} \right) - \arctan\left( \frac{E_\mathrm{F} - \varepsilon - (1/2+a) eV}{\Gamma} \right) \right]. \f]
 */
class SymOneSiteSimulateModel : public SimulateModel,
	public DifferentialConductance, public StaticConductance {

private:
	/**
	 * \brief Ordered list (vector) of the parameters needed for this model.
	 *
	 * If this order is changed, SymOneSiteSimulateModel::unpack_parameters
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
	 * \param[out] gamma The site-lead coupling.
	 * \param[out] a The voltage drop parameter.
	 */
	static void unpack_parameters(const std::vector<double> &vec, double &ef,
		double &v, double &epsilon, double &gamma, double &a);

public:
	SymOneSiteSimulateModel() = delete;
	virtual ~SymOneSiteSimulateModel() = default;

	/**
	 * \internal
	 * \brief Default constructor.
	 *
	 * The constructor specifies the required parameters.
	 *
	 * \param[in] avail The available distributions, keyed by name.
	 * \endinternal
	 */
	SymOneSiteSimulateModel(
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
	 * \param[in] gamma The level-electrode coupling.
	 * \param[in] a The voltage drop scaling factor.
	 * \return The transmission for this set of parameters.
	 */
	static double transmission(const double e, const double v, const double eps,
		const double gamma, const double a);

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
