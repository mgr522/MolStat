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

#ifndef __single_molecule_cv_simulate_model_h__
#define __single_molecule_cv_simulate_model_h__

#include <memory>
#include <string>
#include <vector>
#include <map>
#include <general/random_distributions/rng.h>
#include <general/simulator_tools/simulate_model_interface.h>
#include "transport_observables.h"

#include <sundials/sundials_types.h>
#include <sundials/sundials_dense.h>
#include <cvode/cvode.h>
#include <nvector/nvector_serial.h>
#include <cvode/cvode_dense.h>

#include <gsl/gsl_const_mksa.h>
#include <gsl/gsl_sf_exp.h>
#include <gsl/gsl_math.h>
using std::shared_ptr;

/**
 * \brief Simulator model for single molecule cyclic voltammogram.
 *
 * Model parameters are
 * - `e0` (\f$E_\mathrm{0}\f$), the initial applied potential at \f$t=0\f$,
 * - `eref` (\f$E_\mathrm{Ref}\f$), the reference potential,
 * - `lambda` (\f$\lambda\f$), the reorganization energy,
 * - `af` (\f$A_\mathrm{f}\f$), the prefactor for forward half-reaction rate constant,
 * - `ab` (\f$A_\mathrm{b}\f$), the prefactor for backward half-reaction rate constant,
 * - `v` (\f$v\f$), the sweeping rate of the applied potential,
 * - 'temperature' (\f$T\f$), the temperature,
 * - 'n' (\f$n\f$), the number of electrons involved in the half-reaction,
 * - 'tlimit' (\f$t_\mathrm{lim}\f$), the time instant at which the potential start to decrease.
 *
 * The probabilities \f$P_\mathrm{O}(t)\f$ and f$P_\mathrm{R}(t)\f$ for oxidized species \f$O\f$ and reduced species \f$R\f$ in t:
 * \f$ \hat{\Sigma}_\mathrm{L/R} = \left[ -i\Gamma_\mathrm{L,R}/2 \right] \f$,
 * the transmission function is
 * \f[ T(E) = \frac{4\Gamma_\mathrm{L} \Gamma_\mathrm{R}}{4(E-\varepsilon-aeV)^2 + (\Gamma_\mathrm{L} + \Gamma_\mathrm{R})^2}. \f]
 * - Differential conductance:
 *   \f[ G_\mathrm{d}(V) = \frac{2e^2}{h} \left[ (1/2-a) T(E_\mathrm{F} + eV/2) + (1/2+a) T(E_\mathrm{F} - eV/2) \right]. \f]
 * - Static conductance:
 *   \f[ G_\mathrm{s}(V) = \frac{2e^2}{h} \frac{2\Gamma_\mathrm{L} \Gamma_\mathrm{R}}{eV(\Gamma_\mathrm{L} +\Gamma_\mathrm{R})} \left[ \arctan\left( \frac{2[E_\mathrm{F} - \varepsilon + (1/2-a) eV]}{\Gamma_\mathrm{L} + \Gamma_\mathrm{R}} \right) - \arctan\left( \frac{2[E_\mathrm{F} - \varepsilon - (1/2+a)eV]}{\Gamma_\mathrm{L} + \Gamma_\mathrm{R}} \right) \right]. \f]
 */
class SingleMoleculeCV: public SimulateModel, public SingMolCVPeak {
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
	 * \param[out] e0 The initial applied potential.
	 * \param[out] v The voltage.
	 * \param[out] epsilon The site energy.
	 * \param[out] gammal The left site-lead coupling.
	 * \param[out] gammar The right site-lead coupling.
	 * \param[out] a The voltage drop parameter.
	 */
  static void unpack_parameters(const std::vector<double> &vec, double &e0,
      double &eref, double &lambda, double &af, double &ab, double &v,
      double &n, double &poinitial, double &temperature, double &tlimit, double &test);
//	static void unpack_parameters(const std::vector<double> &vec, double &ef,
//		double &v, double &epsilon, double &gammal, double &gammar, double &a, double &b);

public:
	SingleMoleculeCV() = delete;
	virtual ~SingleMoleculeCV() = default;

	/**
	 * \internal
	 * \brief Constructor specifying the required parameters.
	 *
	 * \param[in] avail The available distributions, keyed by name.
	 * \endinternal
	 */
	SingleMoleculeCV(
		const std::map<std::string, shared_ptr<RandomDistribution>> &avail);

		/**
	 * \brief Returns current peak potentials for a randomly-generated set of model parameters.
	 *    model parameters.
	 * 
	 * \param[in] r The GSL random number generator handle.
	 * \return The static conductance.
	 */
	virtual std::array<double, 2> PeakPotentials(shared_ptr<gsl_rng> r) const
		override;

    /**
     * \brief Calculate the forward half-reaction rate constant for a set of model parameters.
     *
     * \param[in] t time
     * \param[in] vec The vector of model parameters.
     * \return The forward half-reaction rate constant for this set of model parameters.
     */
  static double kf(double t, const std::vector<double> &vec);
    
  static double peak_potentials(const std::vector<double> &vec);

  static double kb(double t, const std::vector<double> &vec);

  static double E_applied(double t, const std::vector<double> &vec);

  static int f(double t, N_Vector y, N_Vector ydot, void *user_data);
  static int g(double t, N_Vector y, double *gout, void *user_data);
  static int Jac(long int N, double t, N_Vector y, N_Vector fy, DlsMat J, void *user_data, 
    N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);
  static int display_parameters(const std::vector<double> &vec);
};

#endif
