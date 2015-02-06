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

#include <string>
#include <vector>
#include <valarray>
#include "observables.h"

// CVODE headers
#include <sundials/sundials_types.h>
#include <sundials/sundials_dense.h>
#include <cvode/cvode.h>
#include <nvector/nvector_serial.h>
#include <cvode/cvode_dense.h>

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
 * - Sweep rate is measured in \f$ k_\mathrm{B} T / (ne)\f$ per second.
 *
 * Model parameters are
 * - `lambda` (\f$\lambda\f$), the reorganization energy,
 * - `af` (\f$A_\mathrm{f}\f$), prefactor for the forward half-reaction rate
 *   constant,
 * - `ab` (\f$A_\mathrm{b}\f$), prefactor for the backward half-reaction rate
 *   constant,
 * - `eref` (\f$E_\mathrm{ref}\f$), the reference potential,
 * - `e0` (\f$E_\mathrm{0}\f$), the initial potential at \f$t=0\f$,
 * - `v` (\f$v\f$), the sweep rate of the applied potential,
 * - `tlim` (\f$t_\mathrm{lim}\f$), the time where the backward sweep begins.
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
 * \f[k_\mathrm{b}(t)=A_\mathrm{b} e^{-\frac{[ne(E(t)-E_\mathrm{ref})-\lambda]^2}{4\lambda k_\mathrm{B}T}}\f]
 * and the potential wave form \f$E(t)\f$ is taken here as \f$E_0+vt\f$ for
 * \f$0\leq t \leq t_\mathrm{lim}\f$ and \f$E_0+2vt_\mathrm{lim}-vt\f$ for
 * \f$t_\mathrm{lim}\leq t\leq 2t_\mathrm{lim}\f$.
 *
 * Additoinally, the expressions of rate constants simplify in the reduced unit system as:
 * \f[k_\mathrm{f}(t)=A_\mathrm{f} e^{-\frac{[E(t)-E_\mathrm{ref}+\lambda]^2}{4\lambda}}\f]
 * \f[k_\mathrm{b}(t)=A_\mathrm{b} e^{-\frac{[E(t)-E_\mathrm{ref}-\lambda]^2}{4\lambda}}\f]

 * The forward electron transfer potential is \f$E(t)\f$ for the
 * \f$ 0 \le t \le t_\mathrm{lim} \f$ satisfying \f$P_\mathrm{O}(t)=1/2\f$
 * (if such a \f$t\f$ exists). The backward potential is similarly defined for
 * \f$ t_\mathrm{lim} < t \le 2 t_\mathrm{lim} \f$.
 */
class NonNernstianReaction :
	public ForwardETPotential,
	public BackwardETPotential
{
public:
	/// Container index for the reorganization energy.
	static const std::size_t Index_lambda;

	/// Container index for the forward half-reaction rate constant.
	static const std::size_t Index_Af;

	/// Container index for the backward half-reaction rate constant.
	static const std::size_t Index_Ab;

	/// Container index for the reference potential.
	static const std::size_t Index_Eref;

	/// Container index for the initial potential.
	static const std::size_t Index_E0;

	/// Container index for the potential sweep speed.
	static const std::size_t Index_v;

	/// Container index for the potential turn-around time.
	static const std::size_t Index_tlim;

protected:
	virtual std::vector<std::string> get_names() const override;

public:
	/// Constructor that initializes the CVODE memory.
	NonNernstianReaction();

	/// Destructor that frees the CVODE memory.
	virtual ~NonNernstianReaction();

	/**
 	 * \brief The forward half-reaction rate constant for a set of model
	 *    parameters.
	 *
	 * \param[in] t The time
	 * \param[in] params The set of model parameters.
	 * \return The forward half-reaction rate constant for this set of model
	 *    parameters.
	 */
	static double kf(double t, const std::valarray<double> &params);

	/**
 	 * \brief The backward half-reaction rate constant for a set of model
	 *    parameters.
	 *
	 * \param[in] t The time
	 * \param[in] params The set of model parameters.
	 * \return The backward half-reaction rate constant for this set of model
	 *    parameters.
	 */
	static double kb(double t, const std::valarray<double> &params);

	/**
	 * \brief The applied potential at the specified time (for a set of model
	 *    parameters).
	 * 
	 * \param[in] t The time.
	 * \param[in] params The set of model parameters.
	 * \return The applied potential.
	 */
	static double E_applied(double t, const std::valarray<double> &params);

	/**
	 * \brief The redox potential where electron transfers during the forward 
	 *    potential sweeping.
	 * 
	 * \param[in] params The set of model parameters.
	 * \return The forward redox potential.
	 */
	virtual double ForwardETP(const std::valarray<double> &params) const
		override;
	/**
	 * \brief The redox potential where electron transfers during the backward 
	 *    potential sweeping.
	 * 
	 * \param[in] params The set of model parameters.
	 * \return The backward redox potential.
	 */
	virtual double BackwardETP(const std::valarray<double> &params) const
		override;

private:
	/// CVODE integration workspace.
	void *cvode_mem;

	/// CVODE vector for storing \f$P_\mathrm{O}(t)\f$ (for some \f$t\f$).
	N_Vector po;

	/**
	 * \brief Function for CVODE that specifies the right-hand side of the
	 *    differential equation for \f$P_\mathrm{O}(t)\f$.
	 * 
	 * \param[in] t The current time.
	 * \param[in] po The current \f$P_\mathrm{O}(t)\f$ (as a `N_Vector`).
	 * \param[out] podot The value \f$P'_\mathrm{O}(t)\f$ (as a `N_Vector`).
	 * \param[in] voidparams The model parameters, in `void*` form (per CVODE).
	 * \return 0 if successful, or a non-zero value if an error occured.
	 */
	static int ode(double t, N_Vector po, N_Vector podot, void *voidparams);

	/**
	 * \brief Function for CVODE that defines the roots to search for; that is,
	 *    times where \f$P_\mathrm{O}(t)=1/2\f$.
	 * 
	 * \param[in] t The current time.
	 * \param[in] po The current \f$P_\mathrm{O}(t)\f$ (as a `N_Vector`).
	 * \param[out] rootvals The value \f$P_\mathrm{O}(t)-1/2\f$, in search of a
	 *    root.
	 * \param[in] voidparams The model parameters, in `void*` form (per CVODE).
	 * \return 0 if successful, or a non-zero value if an error occured.
	 */
	static int half_finder(double t, N_Vector po, double *rootvals,
		void *voidparams);
 
	/**
	 * \brief Function for CVODE that computes the Jacobian, for the ODE.
	 *
	 * \param[in] N The number of differential equations, should be 1 in this
	 *    case.
	 * \param[in] t The current time.
	 * \param[in] po The current \f$P_\mathrm{O}(t)\f$ (as a `N_Vector`).
	 * \param[in] podot The current \f$P_\mathrm{O}'(t)\f$ (as a `N_Vector`).
	 * \param[out] J The dense Jacobian matrix.
	 * \param[in] voidparams The model parameters, in `void*` form (per CVODE).
	 * \param tmp1 A temporary workspace.
	 * \param tmp2 A temporary workspace.
	 * \param tmp3 A temporary workspace.
	 * \return 0 if successful, or a non-zero value if an error occured.
	 */
 	static int jac(long int N, double t, N_Vector po, N_Vector podot, DlsMat J,
 		void *voidparams, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);
};

} // namespace molstat::echem
} // namespace molstat

#endif
