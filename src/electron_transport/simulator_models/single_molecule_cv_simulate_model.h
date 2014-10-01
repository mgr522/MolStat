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
#include <nvector/nvector_serial.h>
#include <sundials/sundials_dense.h>
#include <cvode/cvode.h>
#include <cvode/cvode_dense.h>

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
class SingleMoleculeCV: public SimulateModel, public DifferentialConductance, public SingMolCVPeak {
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
	static void unpack_parameters(const std::vector<double> &vec, double &e0,
		double &eref, double &lambda, double &af, double &ab, double &v, double &n, double &poinitial, double &temperature, double &tlimit);

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
	 * \brief Calculates the differential conductance for a set of model
	 *    parameters.
	 *
	 * \param[in] vec The vector of model parameters.
	 * \return The differential conductance for this set of parameters.
	 */
	static double diff_conductance(const std::vector<double> &vec);
    /**
     * \brief Calculate the forward half-reaction rate for a set of model parameters.
     *
     * \param[in] t The time
     * \param[in] vec The vector of model parameters.
     * \return The forward half-reaction rate for this set of model parameters.
     */
    static double kf(double t, std::vector<double> &vec);
    /**
     * \brief Calculate the backward half-reaction rate for a set of model parameters.
     *
     * \param[in] t The time.
     * \param[in] vec The vector of model parameters.
     * \return The backward half-reaction rate for this set of model parameters.
     */
    static double kb(double t, std::vector<double> &vec);
    /**
     * \brief Calculate the potential applied on the moleculei for a set of model parameters at time t.
     * \param[in] t The time.
     * \param[in] vec The vector of model parameters.
     * \return The applied potential on the molecule at time t.
     */
    static double E_applied(double t, std::vector<double> &vec);
    /**
     * \brief The funtion called by the solver. Calculates the values on the right side of the equations.
     * \param[in] t the time.
     * \param[in] y the vector of functions that is solved.
     */
    static int f(double t, N_Vector y, N_Vector ydot, void *user_data);
    static int g(double t, N_Vector y, double *gout, void *user_data);
    static int Jac(long int N, double t, N_Vector y, N_Vector fy, DlsMat J, void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);
};

#endif
