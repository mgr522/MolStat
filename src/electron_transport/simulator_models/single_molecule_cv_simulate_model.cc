/* This file is a part of MolStat, which is distributed under the Creative
   Commons Attribution-NonCommercial 4.0 International Public License.
   MolStat (c) 2014, Northwestern University. */
/**
 * \file asym_one_site_simulate_model.cc
 * \brief Tight-binding model with one site that couples asymmetrically to
 *    both electrodes.
 *
 * \author Matthew G.\ Reuter
 * \date September 2014
 */

#include "single_molecule_cv_simulate_model.h"
#include <cmath>

#define Ith(v,i)    NV_Ith_S(v,i-1)
#define IJth(A,i,j) DENSE_ELEM(A,i-1,j-1)

#define NEQ 2
#define RTOL 1.0e-5
#define ATOL1 1.0e-8
#define ATOL2 1.0e-8
#define T0 0.0
#define NOUT 1.0

using namespace std;

// if the order of the following list is changed, the unpack_parameters
// function MUST also be updated
const vector<string> SingleMoleculeCV::parameters =
    {"e0", "eref", "lambda", "af", "ab", "v", "n",
    "poinitial", "temperature", "tlimit"};

void SingleMoleculeCV::unpack_parameters(const std::vector<double> &vec,
  double &e0, double &eref, double &lambda, double &af,
  double &ab, double &v, double &n, double &poinitial,
  double &temperature, double &tlimit) {

	e0 = vec[0];
	eref = vec[1];
	lambda = vec[2];
	af = vec[3];
	ab = vec[4];
	v = vec[5];
  n = vec[6];
  poinitial = vec[7];
  temperature = vec[8];
  tlimit = vec[9];
}

SingleMoleculeCV::SingleMoleculeCV(
	const std::map<std::string, shared_ptr<RandomDistribution>> &avail)
	: SimulateModel(avail, parameters) {
}



std::array<double, 2> SingleMoleculeCV::PeakPotentials(shared_ptr<gsl_rng> r)
	const {

	vector<double> params(10);
	double e0, eref, lambda, af, ab, v, n, poinitial, temperature, tlimit;

	// generate and unpack the parameters
	sample(r, params);
	unpack_parameters(params, e0, eref, lambda, af, ab, v, n, poinitial, temperature, tlimit);

	return {{1.0 , peak_potentials(params)}};
}

double SingleMoleculeCV::peak_potentials(const std::vector<double> &vec) {

  double e0, eref, lambda, af, ab, v, n, poinitial, temperature, tlimit;

	// unpack the model parameters
	unpack_parameters(vec, e0, eref, lambda, af, ab, v, n, poinitial, temperature, tlimit);

  double reltol, t, tout;
  N_Vector y, abstol;
  void *cvode_mem;
  int flag, flagr, iout;
  int rootsfound[2];

  std::vector<double> vec_cvode;
  vec_cvode = vec;

  y = abstol = NULL;
  cvode_mem = NULL;

  //Create serial vector of length NEG for I.C and abstol
  y = N_VNew_Serial(NEQ);
  abstol = N_VNew_Serial(NEQ);

  // initialize y
  Ith(y,1) = poinitial;
  Ith(y,2) = 1.0 - poinitial;

  // set the scalar relative tolerance
  reltol = RTOL;

  // set the vector absolute tolerance
  Ith(abstol,1) = ATOL1;
  Ith(abstol,2) = ATOL2;

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

    if (flag == CV_ROOT_RETURN) {
      flagr = CVodeGetRootInfo(cvode_mem, rootsfound);
      root.push_back( E_applied(t, vec) );
    }
    if (flag == CV_SUCCESS) break;
  }
  
  // free y and abstol vectors
  N_VDestroy_Serial(y);
  N_VDestroy_Serial(abstol);

  // free integrator memory
  CVodeFree(&cvode_mem);

  return root[1];
}

double SingleMoleculeCV::kf( double t,
	const std::vector<double> &vec) {

	double e0, eref, lambda, af, ab, v, n, poinitial, temperature, tlimit;

	// unpack the model parameters
	unpack_parameters(vec, e0, eref, lambda, af, ab, v, n, poinitial, temperature, tlimit);

	return af * gsl_sf_exp( - gsl_pow_2( n * GSL_CONST_MKSA_ELECTRON_CHARGE * (E_applied(t, vec) - eref) + lambda * GSL_CONST_MKSA_ELECTRON_CHARGE)
        / (4.0 * lambda * GSL_CONST_MKSA_ELECTRON_CHARGE * GSL_CONST_MKSA_BOLTZMANN * temperature));
}

double SingleMoleculeCV::kb( double t,
	const std::vector<double> &vec) {

	double e0, eref, lambda, af, ab, v, n, poinitial, temperature, tlimit;

	// unpack the model parameters
	unpack_parameters(vec, e0, eref, lambda, af, ab, v, n, poinitial, temperature, tlimit);

	return ab * gsl_sf_exp( - gsl_pow_2( n * GSL_CONST_MKSA_ELECTRON_CHARGE * (E_applied(t, vec) - eref) - lambda * GSL_CONST_MKSA_ELECTRON_CHARGE)
        / (4.0 * lambda * GSL_CONST_MKSA_ELECTRON_CHARGE * GSL_CONST_MKSA_BOLTZMANN * temperature));
}

double SingleMoleculeCV::E_applied(double t,
  const std::vector<double> &vec) {

  double e0, eref, lambda, af, ab, v, n, poinitial, temperature, tlimit;

  //upack the model paramters
  unpack_parameters(vec, e0, eref, lambda, af, ab, v, n, poinitial, temperature, tlimit);

  double E;

  if (t >= 0 && t <= tlimit)
      E = e0 + v * t;
  if (t > tlimit )
      E = e0 + 2.0 * v * tlimit - v * t;
  return E;
}

int SingleMoleculeCV::f(double t, N_Vector y, N_Vector ydot, 
  void *user_data) {

  double y1, y2;
  std::vector<double> *p;
  p = (std::vector<double> *) user_data;

  y1 = Ith(y,1);
  y2 = Ith(y,2);

  Ith(ydot,1) = - kf(t, *p) * y1 + kb(t, *p) *y2;
  Ith(ydot,2) =   kf(t, *p) * y1 - kb(t, *p) *y2;
  
  return(0);
}

int SingleMoleculeCV::g(double t, N_Vector y, double *gout, 
  void *user_data) {

  double y1;

  y1 = Ith(y,1);
  gout [0] = y1 - 0.5;

  return(0);
}

int SingleMoleculeCV::Jac(long int N, double t, N_Vector y,
              N_Vector fy, DlsMat J, void *user_data,
              N_Vector tmp1, N_Vector tmp2, N_Vector tmp3) {
  std::vector<double> *p;

  p = (std::vector<double> *) user_data;

  IJth(J,1,1) = - kf(t, *p);
  IJth(J,1,2) =   kb(t, *p);
  IJth(J,2,1) =   kf(t, *p);
  IJth(J,2,2) = - kb(t, *p);

  return(0);
}



