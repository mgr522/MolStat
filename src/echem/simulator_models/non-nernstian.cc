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

#define Ith(v,i)    NV_Ith_S(v,i-1)
#define IJth(A,i,j) DENSE_ELEM(A,i-1,j-1)

#define NEQ 1
#define RTOL 1.0e-8
#define ATOL1 1.0e-6
#define ATOL2 1.0e-6
#define T0 0.0

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
};

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
  double e0, eref, lambda, af, ab, v, n, poinitial, temperature, tlimit, direction;

  double reltol, t, tout;
  N_Vector y, abstol;
  void *cvode_mem;
  int flag, flagr, iout;
  int rootsfound[2];

  //set the maximum steps for ODE solver
  int mxsteps = 2000;

  //set the maximum step size
  double hmax = 4.0 * tlimit / mxsteps;
  //double hmin = tlimit * 1.0e-12;

  std::vector<double> vec_cvode;
  vec_cvode = vec;

  y = abstol = NULL;
  cvode_mem = NULL;

  //Create serial vector of length NEG for I.C and abstol
  y = N_VNew_Serial(NEQ);
  abstol = N_VNew_Serial(NEQ);

  // initialize y
  Ith(y,1) = poinitial;

  // set the scalar relative tolerance
  reltol = RTOL;

  // set the vector absolute tolerance
  Ith(abstol,1) = ATOL1;

  // call CVodeCreate to create the solver memory and specify the 
  // Backward Differentiation Formula and the use of a Newton iteration.
  cvode_mem = CVodeCreate(CV_BDF, CV_NEWTON);

  // call CVodeInit to initialize the integrator memory and specify the
  // user's right hand side function in y'=f(t,y), the initial time T0, and
  // the initial dependent variable vector y.
  CVodeInit(cvode_mem, f, T0, y);

  // call CVodeSVtolerances to specify the scalar relative tolerance
  // and vector absolute  tolerance
  CVodeSVtolerances(cvode_mem, reltol, abstol);

  // call CVodeSetUserData to pass parameters to user_data;
  CVodeSetUserData(cvode_mem, &vec_cvode);

  // call CVodeSetMaxNumSteps to specify the maximun number of steps to be taken.
  CVodeSetMaxNumSteps(cvode_mem, mxsteps);

  // call CVodeSetMaxStep to specify the maximum step size.
  CVodeSetMaxStep(cvode_mem, hmax);

  // Call CVodeSetMaxErrTestFails to specify the maximun number of error test failures
  // permitted in attempting one step
  CVodeSetMaxErrTestFails(cvode_mem, 20);

  // call CVodeRootInit to specify the root function g with 1 component
  CVodeRootInit(cvode_mem, 1, g);

  // call CVDense to specify the CVDENSE dense linear solver
  CVDense(cvode_mem, NEQ);

  // set the Jacobian routine to Jac
  CVDlsSetDenseJacFn(cvode_mem, Jac);

  // define a vector to store roots.
  std::vector<double> root(1);
  
  tout = 2.0 * tlimit;

  while(1) {

    flag = CVode(cvode_mem, tout, y, &t, CV_NORMAL);
    //PrintOutput(t, Ith(y,1), 1.0 - Ith(y,1));

    if (flag == CV_ROOT_RETURN) {
      flagr = CVodeGetRootInfo(cvode_mem, rootsfound);
      root.push_back( E_applied(t, vec) );
    }
    if (flag == CV_SUCCESS) break;
    //if (flag == CV_TOO_MUCH_WORK) exit(-1);
  }
  
  // free y and abstol vectors
  N_VDestroy_Serial(y);
  N_VDestroy_Serial(abstol);

  // free integrator memory
  CVodeFree(&cvode_mem);

  //printf("the size of the vector is %4d\n", root.size());
  //printf("the first root is %14.8e\n", root[1]);
  //printf("the second root is %14.8e\n", root[2]);

  if (direction == 1.0) return root[1];
  if (direction == 2.0) return root[2];

  return 0;
}

#if 0
int SingleMoleculeCV::f(double t, N_Vector y, N_Vector ydot, 
  void *user_data) {

  double y1;
  std::vector<double> *p;
  p = (std::vector<double> *) user_data;

  y1 = Ith(y,1);

  Ith(ydot,1) = - ( kf(t, *p)  + kb(t, *p) ) * y1 + kb(t, *p);

  //PrintOutput(t, y1,y2);
  
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