#include <stdio.h>
#include <math.h>
#include <memory>
/* header files of GSL */
#include <gsl/gsl_const_mksa.h>
#include <gsl/gsl_sf_exp.h>
#include <gsl/gsl_math.h>

/* header files of CVODE */
#include <cvode/cvode.h>             /* prototypes for CVODE fcts., consts. */
#include <nvector/nvector_serial.h>  /* serial N_Vector types, fcts., macros */
#include <cvode/cvode_dense.h>       /* prototype for CVDense */
#include <sundials/sundials_dense.h> /* definitions DlsMat DENSE_ELEM */
#include <sundials/sundials_types.h> /* definition of type realtype */

#include "PO.h"

/* These macros are defined in order to write code which exactly matches
   the mathematical problem description given above.

   Ith(v,i) references the ith component of the vector v, where i is in
   the range [1..NEQ] and NEQ is defined below. The Ith macro is defined
   using the N_VIth macro in nvector.h. N_VIth numbers the components of
   a vector starting from 0.

   IJth(A,i,j) references the (i,j)th element of the dense matrix A, where
   i and j are in the range [1..NEQ]. The IJth macro is defined using the
   DENSE_ELEM macro in dense.h. DENSE_ELEM numbers rows and columns of a
   dense matrix starting from 0. */

#define Ith(v,i)    NV_Ith_S(v,i-1)       /* Ith numbers components 1..NEQ */
#define IJth(A,i,j) DENSE_ELEM(A,i-1,j-1) /* IJth numbers rows,cols 1..NEQ */


/* Problem Constants */

#define NEQ   2                /* number of equations  */
#define PO_initial    0.0      /* initial y components */
#define PR_initial    1.0
#define RTOL  1.0e-5   /* scalar relative tolerance            */
#define ATOL1 1.0e-8   /* vector absolute tolerance components */
#define ATOL2 1.0e-8
#define T0    0.0      /* initial time           */
#define TMULT 10.0    /* output time factor     */
#define TADD  1.0
#define NOUT  1              /* number of output times */

#define E0_c      0.0     /* initial applied potential */
#define Eref_c    1.130817   /* reference potential */
#define lamb_c    0.6959579 /* reorganization energy */
#define Af_c      5.019018     /* prefactor for forward half-action rate constant */
#define Ab_c      5.034298     /* prefactor for backward half-action rate constant */
#define v_c       0.01    /* Sweeping rate of the applied potential */
#define Temp_c    300.0   /* Temperature */
#define ne_c      1.0     /* Number of electrons involved in the reaction */
#define t1_c      250.0     /* the time from which the potential start to decrease */
#define T1        2.0*t1_c /* first output time      */
using namespace std;

static void unpack_parameters(const std::vector<double> &vec,
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


static double kf( double t,
	const std::vector<double> &vec) {

	double e0, eref, lambda, af, ab, v, n, poinitial, temperature, tlimit;

	// unpack the model parameters
	unpack_parameters(vec, e0, eref, lambda, af, ab, v, n, poinitial, temperature, tlimit);

	return af * gsl_sf_exp( - gsl_pow_2( n * GSL_CONST_MKSA_ELECTRON_CHARGE * (E_applied(t, vec) - eref) + lambda * GSL_CONST_MKSA_ELECTRON_CHARGE )
        / (4.0 * lambda * GSL_CONST_MKSA_ELECTRON_CHARGE * GSL_CONST_MKSA_BOLTZMANN * temperature));
}

static double kb( double t,
	const std::vector<double> &vec) {

	double e0, eref, lambda, af, ab, v, n, poinitial, temperature, tlimit;

	// unpack the model parameters
	unpack_parameters(vec, e0, eref, lambda, af, ab, v, n, poinitial, temperature, tlimit);

	return ab * gsl_sf_exp( - gsl_pow_2( n * GSL_CONST_MKSA_ELECTRON_CHARGE * (E_applied(t, vec) - eref) - lambda * GSL_CONST_MKSA_ELECTRON_CHARGE )
        / (4.0 * lambda * GSL_CONST_MKSA_ELECTRON_CHARGE * GSL_CONST_MKSA_BOLTZMANN * temperature));
}

static double E_applied(double t,
    const std::vector<double> &vec) {

    double e0, eref, lambda, af, ab, v, n, poinitial, temperature, tlimit;

    //upack the model paramters
    unpack_parameters(vec, e0, eref, lambda, af, ab, v, n, poinitial, temperature, tlimit);

    double E;

    if (t >= 0 && t <= tlimit)
        E = e0 + v * t;
    if (t > tlimit && t <= 2.0 * tlimit)
        E = e0 + 2.0 * v * tlimit - v * t;
    return E;
}


int main()
{
    vector<double> vec(10);
    vec[0] = E0_c;
    vec[1] = Eref_c;
    vec[2] = lamb_c;
    vec[3] = Af_c;
    vec[4] = Ab_c;
    vec[5] = v_c;
    vec[6] = ne_c;
    vec[7] = PO_initial;
    vec[8] = Temp_c;
    vec[9] = t1_c;

    printf("kf = %14.6e \n", kf(1.158025e3, vec));
    printf("kb = %14.6e \n", kb(1.158025e3, vec));

    double reltol, t, tout;
    N_Vector y, abstol;
    void *cvode_mem;
    int flag, flagr, iout;
    int rootsfound[2];

    y = abstol = NULL;
    cvode_mem = NULL;

    /* Create serial vector of length NEQ for I.C. and abstol */
    y = N_VNew_Serial(NEQ);
    abstol = N_VNew_Serial(NEQ);

    /* Initialize y */
    Ith(y,1) = PO_initial;
    Ith(y,2) = 1.0 - PO_initial;

    /* Set the scalar relative tolerance */
    reltol = RTOL;
    /* Set the vector absolute tolerance */
    Ith(abstol,1) = ATOL1;
    Ith(abstol,2) = ATOL2;

    /* Call CVodeCreate to create the solver memory and specify the
    * Backward Differentiation Formula and the use of a Newton iteration */
    cvode_mem = CVodeCreate(CV_BDF, CV_NEWTON);

    /* Call CVodeInit to initialize the integrator memory and specify the
    * user's right hand side function in y'=f(t,y), the inital time T0, and
    * the initial dependent variable vector y. */
    flag = CVodeInit(cvode_mem, f, T0, y);

    /* Call CVodeSVtolerances to specify the scalar relative tolerance
    * and vector absolute tolerances */
    flag = CVodeSVtolerances(cvode_mem, reltol, abstol);

    flag = CVodeSetUserData(cvode_mem, &vec);

    /* Call CVodeRootInit to specify the root function g with 2 components */
    flag = CVodeRootInit(cvode_mem, 1, g);

    /* Call CVDense to specify the CVDENSE dense linear solver */
    flag = CVDense(cvode_mem, NEQ);

    /* Set the Jacobian routine to Jac (user-supplied) */
    flag = CVDlsSetDenseJacFn(cvode_mem, Jac);
//    CVodeSetInitStep(cvode_mem, 1.0e-4);

    /* In loop, call CVode, print results, and test for error.
     Break out of loop when NOUT preset output times have been reached.  */
    printf(" \n Single molecule cyclic voltammetry \n\n");

    iout = 0;  tout = T1;
    vector<double> root(1);
    while(1) {
        flag = CVode(cvode_mem, tout, y, &t, CV_NORMAL);
        PrintOutput(t, Ith(y,1), Ith(y,2));

        if (flag == CV_ROOT_RETURN) {
            flagr = CVodeGetRootInfo(cvode_mem, rootsfound);
            PrintRootInfo(rootsfound[0]);
            root.push_back(E_applied(t,vec));
        }

        if (check_flag(&flag, "CVode", 1)) break;
        if (flag == CV_SUCCESS) break;
    }

    /* Print some final statistics */
    PrintFinalStats(cvode_mem);
    printf("v1 = %14.6e  v2 = %14.6e \n", root[1],root[2]);
    printf("kf = %14.6e \n", kf(1.0, vec));
    printf("kb = %14.6e \n", kb(1.0, vec));
    printf("E_applied = %14.6e \n", E_applied(1.0, vec));



    /* Free y and abstol vectors */
    N_VDestroy_Serial(y);
    N_VDestroy_Serial(abstol);

    /* Free integrator memory */
    CVodeFree(&cvode_mem);

    return(0);
}

/*
 *-------------------------------
 * Functions called by the solver
 *-------------------------------
 */

/*
 * f routine. Compute function f(t,y).
 */

static int f(double t, N_Vector y, N_Vector ydot, void *user_data)
{
    double y1, y2;
    vector<double> *p;

    p = (vector<double> *) user_data;

    y1 = Ith(y,1); y2 = Ith(y,2);

    Ith(ydot,1) = - kf(t, *p) * y1 + kb(t, *p) * y2;
    Ith(ydot,2) =   kf(t, *p) * y1 - kb(t, *p) * y2;


    return(0);
}

/*
 * g routine. Compute functions g_i(t,y) for i = 0,1.
 */

static int g(double t, N_Vector y, double *gout, void *user_data)
{
  double y1;

  y1 = Ith(y,1);
  gout [0] = y1 - 0.5;

  return(0);
}

/*
 * Jacobian routine. Compute J(t,y) = df/dy. *
 */

static int Jac(long int N, double t,
               N_Vector y, N_Vector fy, DlsMat J, void *user_data,
               N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
    vector<double> *p;

    p = (vector<double> *) user_data;

    IJth(J,1,1) = - kf(t, *p);
    IJth(J,1,2) =   kb(t, *p);
    IJth(J,2,1) =   kf(t, *p);
    IJth(J,2,2) = - kb(t, *p);

  return(0);
}

/*
 *-------------------------------
 * Private helper functions
 *-------------------------------
 */

static void PrintOutput(double t, double y1, double y2)
{
#if defined(SUNDIALS_EXTENDED_PRECISION)
  printf("At t = %0.4Le PO = %14.6Le  PR = %14.6Le\n", t, y1, y2);
#elif defined(SUNDIALS_DOUBLE_PRECISION)
  printf("At t = %0.4le PO = %14.6le  PR = %14.6le\n", t, y1, y2);
#else
  printf("At t = %0.4e  PO = %14.6e  PR = %14.6e\n", t, y1, y2);
#endif

  return;
}

static void PrintRootInfo(int root_f1)
{
  printf("    rootsfound[] = %3d\n", root_f1);

  return;
}

/*
 * Get and print some final statistics
 */

static void PrintFinalStats(void *cvode_mem)
{
  long int nst, nfe, nsetups, nje, nfeLS, nni, ncfn, netf, nge;
  int flag;

  flag = CVodeGetNumSteps(cvode_mem, &nst);
  check_flag(&flag, "CVodeGetNumSteps", 1);
  flag = CVodeGetNumRhsEvals(cvode_mem, &nfe);
  check_flag(&flag, "CVodeGetNumRhsEvals", 1);
  flag = CVodeGetNumLinSolvSetups(cvode_mem, &nsetups);
  check_flag(&flag, "CVodeGetNumLinSolvSetups", 1);
  flag = CVodeGetNumErrTestFails(cvode_mem, &netf);
  check_flag(&flag, "CVodeGetNumErrTestFails", 1);
  flag = CVodeGetNumNonlinSolvIters(cvode_mem, &nni);
  check_flag(&flag, "CVodeGetNumNonlinSolvIters", 1);
  flag = CVodeGetNumNonlinSolvConvFails(cvode_mem, &ncfn);
  check_flag(&flag, "CVodeGetNumNonlinSolvConvFails", 1);

  flag = CVDlsGetNumJacEvals(cvode_mem, &nje);
  check_flag(&flag, "CVDlsGetNumJacEvals", 1);
  flag = CVDlsGetNumRhsEvals(cvode_mem, &nfeLS);
  check_flag(&flag, "CVDlsGetNumRhsEvals", 1);

  flag = CVodeGetNumGEvals(cvode_mem, &nge);
  check_flag(&flag, "CVodeGetNumGEvals", 1);

  printf("\nFinal Statistics:\n");
  printf("nst = %-6ld nfe  = %-6ld nsetups = %-6ld nfeLS = %-6ld nje = %ld\n",
	 nst, nfe, nsetups, nfeLS, nje);
  printf("nni = %-6ld ncfn = %-6ld netf = %-6ld nge = %ld\n \n",
	 nni, ncfn, netf, nge);
}

/*
 * Check function return value...
 *   opt == 0 means SUNDIALS function allocates memory so check if
 *            returned NULL pointer
 *   opt == 1 means SUNDIALS function returns a flag so check if
 *            flag >= 0
 *   opt == 2 means function allocates memory so check if returned
 *            NULL pointer
 */

static int check_flag(void *flagvalue, char const *funcname, int opt)
{
  int *errflag;

  /* Check if SUNDIALS function returned NULL pointer - no memory allocated */
  if (opt == 0 && flagvalue == NULL) {
    fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed - returned NULL pointer\n\n",
	    funcname);
    return(1); }

  /* Check if flag < 0 */
  else if (opt == 1) {
    errflag = (int *) flagvalue;
    if (*errflag < 0) {
      fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed with flag = %d\n\n",
	      funcname, *errflag);
      return(1); }}

  /* Check if function returned NULL pointer - no memory allocated */
  else if (opt == 2 && flagvalue == NULL) {
    fprintf(stderr, "\nMEMORY_ERROR: %s() failed - returned NULL pointer\n\n",
	    funcname);
    return(1); }

  return(0);
}
