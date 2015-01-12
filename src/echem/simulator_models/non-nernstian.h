/* This file is a part of MolStat, which is distributed under the Creative
   Commons Attribution-NonCommercial 4.0 International Public License.
   MolStat (c) 2014, Northwestern University. */
/**
 * \file non-nernstian.h
 * \brief Simulator model for electron transfer with a non-Nernstian reaction.
 *
 * \author Bo Fu, Matthew G.\ Reuter
 * \date December 2014, January 2015
 */

#ifndef __non_nernstian_reaction_h__
#define __non_nernstian_reaction_h__

#include <memory>
#include <string>
#include <vector>
#include <map>
#include <general/random_distributions/rng.h>
#include <general/simulator_tools/simulate_model_interface.h>
#include "observables.h"

#include <sundials/sundials_types.h>
#include <sundials/sundials_dense.h>
#include <cvode/cvode.h>
#include <nvector/nvector_serial.h>
#include <cvode/cvode_dense.h>

#include <gsl/gsl_const_mksa.h>
#include <gsl/gsl_sf_exp.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_log.h>

namespace molstat {
namespace echem {

/**
 * \brief Simulator model for single molecule electron transfer using a
 *    non-Nernstian reaction with \f$n\f$ electrons transferred.
 *
 * Physical parameters include
 * - \f$T\f$, the temperature,
 * - \f$n\f$, the number of electrons transferred in the reaction.
 * .
 * Note that these physical parameters are not directly used by the model
 * (i.e., they do not need to be specified in the input file), but they are
 * used in the following reduced unit system employed by the model and model
 * parameters.
 * - Time is measured in seconds.
 * - Energy is measured in \f$ k_\mathrm{B} T\f$.
 * - Electric potential is measured in \f$ k_\mathrm{B} T / (ne)\f$.
 *
 * Model parameters are
 * - `eref` (\f$E_\mathrm{ref}\f$), the reference potential,
 * - `lambda` (\f$\lambda\f$), the reorganization energy,
 * - `af` (\f$A_\mathrm{f}\f$), prefactor for the forward half-reaction rate
 *   constant,
 * - `ab` (\f$A_\mathrm{b}\f$), prefactor for the backward half-reaction rate
 *   constant,
 * - `e0` (\f$E_\mathrm{0}\f$), the initial potential at \f$t=0\f$,
 * - `v` (\f$v\f$), the sweep rate of the applied potential,
 * - `tinv` (\f$t_\mathrm{lim}\f$), the time where the backward sweep begins.
 *
 * The probabilities \f$P_\mathrm{O}(t)\f$ and \f$P_\mathrm{R}(t)\f$ that the
 * molecule is in its oxidized (O) or reduced species (R), respectively, in the
 * electrochemical half-reaction
 * \f[ O +ne^-\rightleftharpoons R \f]
 * evolve with time by the following master equation
 * \f[P_\mathrm{O}'(t)=-k_\mathrm{f}(t)P_\mathrm{O}(t) + k_\mathrm{b}(t)P_\mathrm{R}(t),\f]
 * \f[P_\mathrm{R}'(t)= k_\mathrm{f}(t)P_\mathrm{O}(t) - k_\mathrm{b}(t)P_\mathrm{R}(t),\f]
 * \f[P_\mathrm{O}+P_\mathrm{R}=1.\f]
 * The initial conditions may be \f$P_\mathrm{O}=1\f$ and \f$P_\mathrm{R}=0\f$.
 *
 * The rate constants of the forward and backward half-reactions are described
 * using Marcus theory,
 * \f[k_\mathrm{f}(t)=A_\mathrm{f} e^{-\frac{[ne(E(t)-E_\mathrm{ref})+\lambda]^2}{4\lambda k_\mathrm{B}T}}\f]
 * \f[k_\mathrm{f}(t)=A_\mathrm{f} e^{-\frac{[ne(E(t)-E_\mathrm{ref})+\lambda]^2}{4\lambda k_\mathrm{B}T}}\f]
 * and the potential wave form \f$E(t)\f$ is taken here as \f$E_0+vt\f$ for
 * \f$0\leq t \leq t_\mathrm{lim}\f$ and \f$E_0+2vt_\mathrm{lim}-vt\f$ for
 * \f$t_\mathrm{lim}\leq t\leq 2t_\mathrm{lim}\f$.
 *
 * The forward electron transfer potential is \f$E(t)\f$ for the
 * \f$ 0 \le t \le t_\mathrm{lim} \f$ satisfying \f$P_\mathrm{O}(t)=1/2\f$
 * (if such a \f$t\f$ exists). The backward potential is similarly defined for
 * \f$ t_\mathrm{lim} < t \le 2 t_\mathrm{lim} \f$.
 */
class SingleMoleculeCV : public SimulateModel, public SingMolCVPeak {
private:
	/**
	 * \brief Ordered list (vector) of the parameters needed for this model.
	 *
	 * If this order is changed, SingleMoleculeCV::unpack_parameters
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
   * \param[in] t The time
   * \param[in] vec The vector of model parameters.
   * \return The forward half-reaction rate constant for this set of model parameters.
   */
  static double kf(double t, const std::vector<double> &vec);

  /**
   * \brief Return the backward half-reaction rate constant for a set of model parameters.
   *
   * \param[in] t The time
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

  static int PrintOutput(double t, double y1, double y2);
};

} // namespace molstat::echem
} // namespace molstat

#endif
