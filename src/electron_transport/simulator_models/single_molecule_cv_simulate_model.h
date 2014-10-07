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
#include <gsl/gsl_sf_log.h>
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
 * - `temperature` (\f$T\f$), the temperature,
 * - `n` (\f$n\f$), the number of electrons involved in the half-reaction,
 * - `tlimit` (\f$t_\mathrm{lim}\f$), the time when the applied potential inverts.
 *
 * The probabilities \f$P_\mathrm{O}(t)\f$ and \f$P_\mathrm{R}(t)\f$ for oxidized species \f$O\f$ and reduced species \f$R\f$ in the
 * electrochemical half-reaction
 * \f[ O +ne^-\rightleftharpoons R \f]
 * evolve with time by the following master equations
 * \f[P_\mathrm{O}'(t)=-k_\mathrm{f}(t)P_\mathrm{O}(t) + k_\mathrm{b}(t)P_\mathrm{R}(t)\f]
 * \f[P_\mathrm{R}'(t)= k_\mathrm{f}(t)P_\mathrm{O}(t) - k_\mathrm{b}(t)P_\mathrm{R}(t)\f]
 * \f[P_\mathrm{O}+P_\mathrm{R}=1\f]
 * with initial conditions \f$P_\mathrm{O}=1.0\f$ or \f$P_\mathrm{O}=0\f$, the rate constants of forward and backward half-reactions described using Marcus theory
 * \f[k_\mathrm{f}(t)=A_\mathrm{f} e^{-\frac{[ne(E(t)-E_\mathrm{Ref})+\lambda]^2}{4\lambda k_\mathrm{B}T}}\f]
 * \f[k_\mathrm{f}(t)=A_\mathrm{f} e^{-\frac{[ne(E(t)-E_\mathrm{Ref})+\lambda]^2}{4\lambda k_\mathrm{B}T}}\f]
 * and the potential wave form \f$E(t)\f$ which is \f$E_0+vt\f$ for \f$0\leq t \leq t_\mathrm{lim}\f$ and \f$E_0+2vt_\mathrm{lim}-vt\f$ for \f$t_\mathrm{lim}\leq t\leq 2t_\mathrm{lim}\f$.
 * - Peak potentials:
 *   The current in single molecule electrochemistry is non-zero only at \f$t_\mathrm{P}\f$ when the half-reaction takes place. We assume that \f$t_\mathrm{P}\f$ are the roots of the 
 *   following equations:
 *   \f[P_\mathrm{O}(t)=0.5\f] 
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
	 * \param[out] eref The reference potential(V).
	 * \param[out] lambda The reorganization energy(eV).
	 * \param[out] af The prefactor for forward half-reaction rate constant.
   * \param[out] ab The prefactor for backward half-reaction rate constant.
	 * \param[out] v The sweeping rate of the applied potential.
	 * \param[out] temperature The temperature of the environment.
   * \param[out] n The number of electrons involved in the half-reaction.
   * \param[out] tlimit The time when the applied potential inverts.
   * \param[out] direction The switch that determine which peak will be returned.
	 */
  static void unpack_parameters(const std::vector<double> &vec, double &e0,
      double &eref, double &lambda, double &af, double &ab, double &v,
      double &n, double &poinitial, double &temperature, double &tlimit, double &direction);
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
     * \brief Return the forward half-reaction rate constant for a set of model parameters.
     *
     * \param[in] t time
     * \param[in] vec The vector of model parameters.
     * \return The forward half-reaction rate constant for this set of model parameters.
     */
  static double kf(double t, const std::vector<double> &vec);
 
    /**
     * \brief Return the backward half-reaction rate constant for a set of model parameters.
     *
     * \param[in] t time
     * \param[in] vec The vector of model parameters.
     * \return The forward half-reaction rate constant for this set of model parameters.
     */
  static double kb(double t, const std::vector<double> &vec);

	/**
	 * \brief Returns current peak potentials for a set of model parameters.
	 * 
	 * \param[in] vec The vector of model parameters.
	 * \return The peak potentials.
	 */
  static double peak_potentials(const std::vector<double> &vec);

	/**
	 * \brief Returns the applied potential for a set of model parameters at time \f$t\f$.
	 * 
     * \param[in] t The time.
	 * \param[in] vec The vecotr of model parameter.
	 * \return The applied potential.
	 */
  static double E_applied(double t, const std::vector<double> &vec);

	/**
	 * \brief Callable function by CVODE. Defines the right side of the differential equations \f$y'=f(t,y)\f$
     * where \f$y=(P_\mathrm{O}(t),P_\mathrm{R}(t)\f$.
	 * 
     * \param[in] t The current value of time.
	 * \param[in] y The current value of vector \f$[P_\mathrm{O}(t),P_\mathrm{R}(t)]\f$.
     * \param[out] ydot The output vector \f$[P'_\mathrm{O}(t),P'_\mathrm{R}(t)\]f$.
	 * \param[in] user_data The pointer passed to CVodeSetUserData.
     * \return 0 if successful or a non-zero value if an error occured.
	 */
  static int f(double t, N_Vector y, N_Vector ydot, void *user_data);

    /**
	 * \brief Callable function by CVODE. Defines the roots to search for.
	 * 
     * \param[in] t The current value of time.
	 * \param[in] y The current value of vector \f$[P_\mathrm{O}(t),P_\mathrm{R}(t)]\f$.
     * \param[out] gout The output vector of roots.
	 * \param[in] user_data The pointer passed to CVodeSetUserData.
     * \return 0 if successful or a non-zero value if an error occured.
	 */
  static int g(double t, N_Vector y, double *gout, void *user_data);
 
    /**
	 * \brief Callable function by CVODE. Computes the dense Jacobian \f$J=\partial f/\partial y\f$.
	 *
     * \param[in] N The number of differential equations.
     * \param[in] t The current value of time.
	 * \param[in] y The current value of vector \f$[P_\mathrm{O}(t),P_\mathrm{R}(t)]\f$.
     * \param[in] fy The current value of vector \f$f(t,y)\f$.
     * \param[out] Jac The dense Jacobian matrix.
	 * \param[in] user_data The pointer passed to CVodeSetUserData.
     * \return 0 if successful or a non-zero value if an error occured.
	 */
  static int Jac(long int N, double t, N_Vector y, N_Vector fy, DlsMat J, void *user_data, 
    N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);
 
    /**
	 * \brief Display all the model parameters in the current run.
	 *
     * \param[in] vec The vector of model parameters. 
     * \return 0 
	 */

  static int display_parameters(const std::vector<double> &vec);
};

#endif
