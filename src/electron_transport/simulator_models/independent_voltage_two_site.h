/* This file is a part of MolStat, which is distributed under the Creative
   Commons Attribution-NonCommercial 4.0 International Public License.
   MolStat (c) 2014, Northwestern University. */
/**
 * \file independent_voltage_one_site.h
 * \brief A sum of two independent, symmetric-coupling, voltage-dependent
 *    one-site tight-binding models for calculating conductances.
 *
 * \author Matthew G.\ Reuter
 * \date July 2014
 */

#ifndef __independent_voltage_two_site_h__
#define __independent_voltage_two_site_h__

#include <memory>
#include <general/random_distributions/rng.h>
#include "conductance_model.h"
#include "symmetric_voltage_one_site.h"

using std::shared_ptr;

/**
 * \brief Class encapsulating transport through two independent channels, where
 *    each channel is a voltage-dependent, one-site model (symmetric coupling).
 *
 * Model parameters are
 * - `epsilon1` (\f$\varepsilon_1\f$), the site-energy for channel 1.
 * - `gamma1` (\f$\Gamma_1\f$), the site/lead coupling for channel 1.
 * - `a1` (\f$a_1\f$), the strength of the voltage dependence for channel 1.
 * - `epsilon2` (\f$\varepsilon_2\f$), the site-energy for channel 2.
 * - `gamma2` (\f$\Gamma_2\f$), the site/lead coupling for channel 2.
 * - `a2` (\f$a_2\f$), the strength of the voltage dependence for channel 2.
 *
 * This model assumes the system has two independent channels, each of which is
 * described by the SymmetricVoltageOneSiteModel. The six parameters above are
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
 * SymmetricVoltageOneSiteModel::transmission) and \f$T_2(E)\f$ is the
 * transmission from channel 2. Accordingly, the differential and static
 * conductances are just the sums of the respective quantities through each
 * channel.
 */
class IndependentVoltageTwoSiteModel : public ConductanceModel {
protected:
	/**
	 * \internal
	 * \brief Model for channel 1.
	 * \endinternal
	 */
	const SymmetricVoltageOneSiteModel channel1;

	/**
	 * \internal
	 * \brief Model for channel 2.
	 * \endinternal
	 */
	const SymmetricVoltageOneSiteModel channel2;

public:
	IndependentVoltageTwoSiteModel() = delete;

	/**
	 * \brief Constructor requiring the random distributions.
	 *
	 * \param[in] eps1 The distribution for \f$\varepsilon_1\f$.
	 * \param[in] gamma1 The distribution for \f$\Gamma_1\f$.
	 * \param[in] a1 The distribution for \f$a_1\f$.
	 * \param[in] eps2 The distribution for \f$\varepsilon_2\f$.
	 * \param[in] gamma2 The distribution for \f$\Gamma_2\f$.
	 * \param[in] a2 The distribution for \f$a_2\f$.
	 */
	IndependentVoltageTwoSiteModel(
		const shared_ptr<const RandomDistribution> &eps1,
		const shared_ptr<const RandomDistribution> &gamma1,
		const shared_ptr<const RandomDistribution> &a1,
		const shared_ptr<const RandomDistribution> &eps2,
		const shared_ptr<const RandomDistribution> &gamma2,
		const shared_ptr<const RandomDistribution> &a2);

	/**
	 * \internal
	 * \brief Destructor.
	 * \endinternal
	 */
	virtual ~IndependentVoltageTwoSiteModel() = default;

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
};

#endif
