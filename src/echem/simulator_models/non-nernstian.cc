/* This file is a part of MolStat, which is distributed under the Creative
   Commons Attribution-NonCommercial 4.0 International Public License.
   MolStat (c) 2014, Northwestern University. */
/**
 * \file non-nernstian.cc
 * \brief Simulator model for electron transfer with a non-Nernstian reaction.
 *
 * \author Bo Fu, Matthew G.\ Reuter
 * \date December 2014, January 2015
 */

#include "non-nernstian.h"
#include <cmath>
#include <general/simulator_tools/simulator_exceptions.h>

/// \cond
#define Ith(v,i)    NV_Ith_S(v,i-1)
#define IJth(A,i,j) DENSE_ELEM(A,i-1,j-1)

#define RTOL 1.e-8
#define ATOL 1.e-6
#define MAXSTEPS 2000
/// \endcond

namespace molstat {
namespace echem {

const std::size_t NonNernstianReaction::Index_lambda = 0;
const std::size_t NonNernstianReaction::Index_Af = 1;
const std::size_t NonNernstianReaction::Index_Ab = 2;
const std::size_t NonNernstianReaction::Index_Eref = 3;
const std::size_t NonNernstianReaction::Index_E0 = 4;
const std::size_t NonNernstianReaction::Index_v = 5;
const std::size_t NonNernstianReaction::Index_tlim = 6;

std::vector<std::string> NonNernstianReaction::get_names() const
{
	std::vector<std::string> ret(7);

	ret[Index_lambda] = "lambda";
	ret[Index_Af] = "af";
	ret[Index_Ab] = "ab";
	ret[Index_Eref] = "eref";
	ret[Index_E0] = "e0";
	ret[Index_v] = "v";
	ret[Index_tlim] = "tlim";

	return ret;
}

double NonNernstianReaction::kf(double t, const std::valarray<double> &params)
{
	// unpack the parameters
	const double &eref = params[Index_Eref];
	const double &lambda = params[Index_lambda];
	const double &af = params[Index_Af];
	
	double exponential = E_applied(t, params) - eref + lambda;
	exponential *= -0.25 * exponential / lambda;

	double log_kf = exponential + log(af);
	if(log_kf < -650.)
		return 0.;
	else
		return exp(log_kf);
}

double NonNernstianReaction::kb(double t, const std::valarray<double> &params)
{
	// unpack the parameters
	const double &eref = params[Index_Eref];
	const double &lambda = params[Index_lambda];
	const double &ab = params[Index_Ab];
	
	double exponential = E_applied(t, params) - eref - lambda;
	exponential *= -0.25 * exponential / lambda;

	double log_kb = exponential + log(ab);
	if(log_kb < -650.)
		return 0.;
	else
		return exp(log_kb);
}

double NonNernstianReaction::E_applied(double t,
	const std::valarray<double> &params)
{
	// unpack the parameters
	const double &e0 = params[Index_E0];
	const double &v = params[Index_v];
	const double &tlim = params[Index_tlim];

	double E;

	if(t >= 0. && t <= tlim)
      E = e0 + v * t;
	else if(t > tlim && t <= 2.*tlim)
      E = e0 + 2.0 * v * tlim - v * t;
	else
      throw molstat::NoObservableProduced();

   return E;
}

NonNernstianReaction::NonNernstianReaction()
	: cvode_mem{nullptr}, po{}
{
	// create serial vectors for the initial condition and for abstol
	po = N_VNew_Serial(1);

	// initialize po (initial condition)
	Ith(po, 1) = 0.;

	// create the solver memory and specify the Backward Differentiation Formula
	// and the use of a Newton iteration.
	cvode_mem = CVodeCreate(CV_BDF, CV_NEWTON);

	// initialize the integrator memory and specify the right hand side function
	// in po'=f(t, po), the initial time 0, and the initial dependent variable
	// vector, po.
	CVodeInit(cvode_mem, ode, 0., po);

	// specify the scalar relative tolerance and vector absolute tolerance
	CVodeSStolerances(cvode_mem, RTOL, ATOL);

	// specify the maximum number of steps to take
	CVodeSetMaxNumSteps(cvode_mem, MAXSTEPS);

	// specify the maximum number of error test failures permitted in attempting
	// one step
	CVodeSetMaxErrTestFails(cvode_mem, 20);

	// specify the CVDENSE dense linear solver
	CVDense(cvode_mem, 1);

	// set the Jacobian routine
	CVDlsSetDenseJacFn(cvode_mem, jac);
}

NonNernstianReaction::~NonNernstianReaction()
{
	N_VDestroy_Serial(po);
	CVodeFree(&cvode_mem);
}

double NonNernstianReaction::ForwardETP(const std::valarray<double> &params)
	const
{
	// unpack the parameters
	const double &tlim = params[Index_tlim]; // sweep turnaround time
	double t, tmax;

	// set the maximum step size
	CVodeSetMaxStep(cvode_mem, 4. * tlim / MAXSTEPS);

	// pass model parameters to the CVODE functions
	std::valarray<double> params_copy { params };
	CVodeSetUserData(cvode_mem, &params_copy);	

	// forward sweep
	// reinitialize the integrator
	Ith(po, 1) = 0.;
	CVodeReInit(cvode_mem, 0., po);

	// set the maximum time for propagation
	tmax = tlim;

	// specify the root function, where we search for the time when po=0.5
	CVodeRootInit(cvode_mem, 1, half_finder);

	int flag;
	do
	{
		flag = CVode(cvode_mem, tmax, po, &t, CV_NORMAL);

		if(flag == CV_ROOT_RETURN)
		{
			// found the time!

			// clean up the integrator
			CVodeSetUserData(cvode_mem, nullptr);
			CVodeRootInit(cvode_mem, 0, nullptr);

			// return the potential
			return E_applied(t, params);
		}
	} while(flag != CV_SUCCESS);

	// we made it to the end of the forward sweep and didn't find a time
	// when po = 0.5 -- no observable to report

	// clean up the integrator
	CVodeSetUserData(cvode_mem, nullptr);
	CVodeRootInit(cvode_mem, 0, nullptr);

	throw molstat::NoObservableProduced();

	// we don't care about the backward sweep
}

double NonNernstianReaction::BackwardETP(const std::valarray<double> &params)
	const
{
	// unpack the parameters
	const double &tlim = params[Index_tlim]; // sweep turnaround time
	double t, tmax;

	// set the maximum step size
	CVodeSetMaxStep(cvode_mem, 4. * tlim / MAXSTEPS);

	// pass model parameters to the CVODE functions
	std::valarray<double> params_copy { params };
	CVodeSetUserData(cvode_mem, &params_copy);	

	// forward sweep
	// reinitialize the integrator
	Ith(po, 1) = 0.;
	CVodeReInit(cvode_mem, 0., po);

	// set the maximum time for propagation
	tmax = tlim;

	int flag;
	do
	{
		flag = CVode(cvode_mem, tmax, po, &t, CV_NORMAL);
	} while(flag != CV_SUCCESS);

	// now do the backward sweep

	// specify the root function, where we search for the time when po=0.5
	CVodeRootInit(cvode_mem, 1, half_finder);

	// set the maximum time for propagation
	tmax = 2.*tlim;

	do
	{
		flag = CVode(cvode_mem, tmax, po, &t, CV_NORMAL);

		if(flag == CV_ROOT_RETURN)
		{
			// found the time!

			// clean up the integrator
			CVodeSetUserData(cvode_mem, nullptr);
			CVodeRootInit(cvode_mem, 0, nullptr);

			// return the potential
			return E_applied(t, params);
		}
	} while(flag != CV_SUCCESS);

	// we made it to the end of the backward sweep and didn't find a time
	// when po = 0.5 -- no observable to report

	// clean up the integrator
	CVodeSetUserData(cvode_mem, nullptr);
	CVodeRootInit(cvode_mem, 0, nullptr);

	throw molstat::NoObservableProduced();
}

int NonNernstianReaction::ode(double t, N_Vector po, N_Vector podot, 
	void *voidparams)
{
	const std::valarray<double> &params =
		*(const std::valarray<double> *)voidparams;

	const double poval = Ith(po, 1);
	const double kfval = kf(t, params);
	const double kbval = kb(t, params);

	Ith(podot, 1) = kbval - (kfval  + kbval) * poval;

	return 0;
}

int NonNernstianReaction::half_finder(double t, N_Vector po, double *rootvals, 
	void *voidparams)
{
	const double poval = Ith(po, 1);
	rootvals[0] = poval - 0.5;

	return 0;
}

int NonNernstianReaction::jac(long int N, double t, N_Vector po,
	N_Vector podot, DlsMat J, void *voidparams, N_Vector tmp1, N_Vector tmp2,
	N_Vector tmp3)
{
	const  std::valarray<double> &params =
		*(const std::valarray<double> *)voidparams;

	IJth(J, 1, 1) = - kf(t, params) - kb(t, params);

	return 0;
}

} // namespace molstat::echem
} // namespace molstat
