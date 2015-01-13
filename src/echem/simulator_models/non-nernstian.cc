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
#include <memory>
#include <general/simulator_tools/simulator_exceptions.h>

// CVODE headers
#include <sundials/sundials_types.h>
#include <sundials/sundials_dense.h>
#include <cvode/cvode.h>
#include <nvector/nvector_serial.h>
#include <cvode/cvode_dense.h>

#define Ith(v,i)    NV_Ith_S(v,i-1)
#define IJth(A,i,j) DENSE_ELEM(A,i-1,j-1)

#define RTOL 1.e-8
#define ATOL 1.e-6
#define MAXSTEPS 2000

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

	double log_kf = e + log(af);
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

	double log_kb = e + log(ab);
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
      E = e0 + 2.0 * v * tlimit - v * t;
	else
      throw molstat::NoObservableProduced();

   return E;
}

double NernstianReaction::ForwardETP(const std::valarray<double> &params)
{
	// unpack the parameters
	const double &tlim = params[Index_tlim];

	double t, tmax;
	std::unique_ptr<N_Vector, typeof(&N_VDestroy_Serial)>
		po { nullptr, &N_VDestroy_Serial },
		abstol { nullptr, &N_VDestroy_Serial };
	void *cvode_mem { nullptr };

	// set the maximum step size
	const double dtmax = 4. * tlim / MAXSTEPS;

	std::valarray<double> params_copy { params };

	// Create serial vectors for the initial condition and for abstol
	po = N_VNew_Serial(1);
	abstol = N_VNew_Serial(1);

	// initialize y (initial condition)
	Ith(y, 1) = 1.;

	// set the scalar relative tolerance
	reltol = RTOL;

	// set the vector absolute tolerance
	Ith(abstol, 1) = ATOL;

	// create the solver memory and specify the Backward Differentiation Formula
	// and the use of a Newton iteration.
	cvode_mem = CVodeCreate(CV_BDF, CV_NEWTON);

	// initialize the integrator memory and specify the right hand side function
	// in po'=f(t, po), the initial time T0, and the initial dependent variable
	// vector, po.
	CVodeInit(cvode_mem, f, 0., po);

	// specify the scalar relative tolerance and vector absolute tolerance
	CVodeSVtolerances(cvode_mem, RTOL, ATOL);

	// pass model parameters to the CVODE functions
	CVodeSetUserData(cvode_mem, &params_copy);

	// specify the maximum number of steps to take
	CVodeSetMaxNumSteps(cvode_mem, MAXSTEPS);

	// specify the maximum step size
	CVodeSetMaxStep(cvode_mem, dtmax);

	// specify the maximum number of error test failures permitted in attempting
	// one step
	CVodeSetMaxErrTestFails(cvode_mem, 20);

	// specify the root function g with 1 component
	CVodeRootInit(cvode_mem, 1, g);

	// specify the CVDENSE dense linear solver
	CVDense(cvode_mem, 1);

	// set the Jacobian routine to Jac
	CVDlsSetDenseJacFn(cvode_mem, Jac);

	// set the maximum time for propagation
	tmax = 2. * tlim;

	while(true)
	{
		int flag = CVode(cvode_mem, tmax, po, &t, CV_NORMAL);

		if(flag == CV_ROOT_RETURN)
		{
			int flagr = CVodeGetRootInfo(cvode_mem, rootsfound);
			// free the integrator memory
			CVodeFree(&cvode_mem);
			return E_applied(t, params);
		}
		else if(flag == CV_SUCCESS)
		{
			// made it all the way to the maximum time and didn't find the
			// potential
			CVodeFree(&cvode_mem);
			throw molstat::NoObservableProduced();
		}
	}

	// should not ever be here
	throw molstat::NoObservableProduced();
	return 0.;
}

#if 0
int SingleMoleculeCV::f(double t, N_Vector y, N_Vector ydot, 
  void *user_data) {

  double y1;
  std::vector<double> *p;
  p = (std::vector<double> *) user_data;

  y1 = Ith(y,1);

  Ith(ydot,1) = - ( kf(t, *p)  + kb(t, *p) ) * y1 + kb(t, *p);
  
  return 0;
}

int SingleMoleculeCV::g(double t, N_Vector y, double *gout, 
  void *user_data) {

  double y1;

  y1 = Ith(y,1);
  gout [0] = y1 - 0.5;

  return 0;
}

int SingleMoleculeCV::Jac(long int N, double t, N_Vector y,
              N_Vector fy, DlsMat J, void *user_data,
              N_Vector tmp1, N_Vector tmp2, N_Vector tmp3) {
  std::vector<double> *p;

  p = (std::vector<double> *) user_data;

  IJth(J,1,1) = - kf(t, *p) - kb(t, *p);

  return 0;
}

#endif

} // namespace molstat::echem
} // namespace molstat