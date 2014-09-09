/* This file is a part of MolStat, which is distributed under the Creative
   Commons Attribution-NonCommercial 4.0 International Public License.
   MolStat (c) 2014, Northwestern University. */
/**
 * \file ind_two_chan_simulate_model.h
 * \brief Generic model for two independent channels that each couple
 *    symmetrically to both electrodes.
 *
 * \author Matthew G.\ Reuter
 * \date September 2014
 */

#ifndef __ind_two_chan_simulate_model_h__
#define __ind_two_chan_simulate_model_h__

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
 * - `epsilon1` (\f$\varepsilon_1\f$), the site-energy for channel 1.
 * - `gamma1` (\f$\Gamma_1\f$), the site/lead coupling for channel 1.
 * - `a1` (\f$a_1\f$), the strength of the voltage dependence for channel 1.
 * - `epsilon2` (\f$\varepsilon_2\f$), the site-energy for channel 2.
 * - `gamma2` (\f$\Gamma_2\f$), the site/lead coupling for channel 2.
 * - `a2` (\f$a_2\f$), the strength of the voltage dependence for channel 2.
 *
 * This model assumes the system has two independent channels, each of which is
 * described by the SymOneSiteSimulateModel. The parameters above are
 * the parameters for each channel within the one-site model.
 *
 * Note that these independent channels are, essentially, what would be
 * obtained from performing an eigenchannel decomposition of the system
 * \cite paulsson-115117.
 *
 * In any event, this model does not necessarily correspond to a specific
 * tight-binding model that leads to a transmission function. Instead,
 * the transmission function is
 * \f[ T(E) = T_1(E) + T_2(E), \f]
 * where \f$T_1(E)\f$ is the transmission from channel 1 (using
 * SymOneSiteSimulateModel::transmission) and \f$T_2(E)\f$ is the
 * transmission from channel 2. Accordingly, the differential and static
 * conductances are just the sums of the respective quantities through each
 * channel.
 */
class IndTwoChanSimulateModel : public SimulateModel,
	public DifferentialConductance, public StaticConductance {

private:
	/**
	 * \brief Ordered list (vector) of the parameters needed for this model.
	 *
	 * If this order is changed, IndTwoChanSimulateModel::unpack_parameters
	 * also needs to be updated.
	 */
	static const std::vector<std::string> parameters;

	/**
	 * \brief Unpack a set of parameters from a vector to doubles.
	 *
	 * \param[in] vec The vector containing a set of parameters.
	 * \param[out] ef The Fermi energy.
	 * \param[out] v The voltage.
	 * \param[out] epsilon1 The site energy for channel 1.
	 * \param[out] gamma1 The site-lead coupling for channel 1.
	 * \param[out] a1 The voltage drop parameter for channel 1.
	 * \param[out] epsilon2 The site energy for channel 2.
	 * \param[out] gamma2 The site-lead coupling for channel 2.
	 * \param[out] a2 The voltage drop parameter for channel 2.
	 */
	static void unpack_parameters(const std::vector<double> &vec, double &ef,
		double &v, double &epsilon1, double &gamma1, double &a1,
		double &epsilon2, double &gamma2, double &a2);

public:
	IndTwoChanSimulateModel() = delete;
	virtual ~IndTwoChanSimulateModel() = default;

	/**
	 * \internal
	 * \brief Default constructor.
	 *
	 * The constructor specifies the required parameters.
	 *
	 * \param[in] avail The available distributions, keyed by name.
	 * \endinternal
	 */
	IndTwoChanSimulateModel(
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
	 * \param[in] eps1 The level energy for channel 1.
	 * \param[in] gamma1 The level-electrode coupling for channel 1.
	 * \param[in] a1 The voltage drop scaling factor for channel 1.
	 * \param[in] eps2 The level energy for channel 2.
	 * \param[in] gamma2 The level-electrode coupling for channel 2.
	 * \param[in] a2 The voltage drop scaling factor for channel 2.
	 * \return The transmission for this set of parameters.
	 */
	static double transmission(const double e, const double v,
		const double eps1, const double gamma1, const double a1,
		const double eps2, const double gamma2, const double a2);

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
