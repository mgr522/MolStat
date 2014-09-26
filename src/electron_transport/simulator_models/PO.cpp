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
#define PO_initial    RCONST(0.0)      /* initial y components */
#define PR_initial    RCONST(1.0)
#define RTOL  RCONST(1.0e-5)   /* scalar relative tolerance            */
#define ATOL1 RCONST(1.0e-8)   /* vector absolute tolerance components */
#define ATOL2 RCONST(1.0e-8)
#define T0    RCONST(0.0)      /* initial time           */
#define TMULT RCONST(10.0)     /* output time factor     */
#define TADD  RCONST(1.0)
#define NOUT  1              /* number of output times */

#define E0_c      RCONST(0.0)     /* initial applied potential */
#define Eref_c    RCONST(1.1)     /* reference potential */
#define lamb_c    RCONST(1.0e-19) /* reorganization energy */
#define Af_c      RCONST(5.0)     /* prefactor for forward half-action rate constant */
#define Ab_c      RCONST(5.0)     /* prefactor for backward half-action rate constant */
#define v_c       RCONST(0.01)    /* Sweeping rate of the applied potential */
#define Temp_c    RCONST(300.0)   /* Temperature */
#define ne_c      RCONST(1.0)     /* Number of electrons involved in the reaction */
#define t1_c      RCONST(250.0)     /* the time from which the potential start to decrease */
#define T1        2.0*t1_c /* first output time      */
using namespace std;

struct function_params params = {E0_c, Eref_c, lamb_c, Af_c, Ab_c, v_c, ne_c, PO_initial, Temp_c, t1_c};

double E_applied (double t, function_params* params)
{
    function_params *p = params;

    double E;

    if (t >= 0 && t <= p->t1)
        E = p->E0 + p->v * t;
    if (t > p->t1 && t <= RCONST(2.0) * p->t1)
        E = p->E0 + RCONST(2.0) * p->v * p->t1 - p->v * t;
    if (t > RCONST(2.0) * p->t1)
    {
        printf("time t is out of range!\n");
    }
    return E;
}

double kf (double t, function_params  *params)
{
    function_params *p = params;

//    realtype E; /* applied potential */
//
//    if (t >= 0 && t <= p->t1)
//        E = p->E0 + p->v * t;
//    if (t > p->t1 && t <= RCONST(2.0) * p->t1)
//        E = p->E0 + RCONST(2.0) * p->v * p->t1 - p->v * t;
//    if (t > RCONST(2.0) * p->t1)
//    {
//        printf("time t is out of range!\n");
//    }

    double result;

    result = p->Af * gsl_sf_exp( - gsl_pow_2( p->n * GSL_CONST_MKSA_ELECTRON_CHARGE * (E_applied(t, params) - p->E_ref) + p->lamb) / (4.0 * p->lamb * GSL_CONST_MKSA_BOLTZMANN * p->T));
    return result;

}
double kb (double t, function_params *params)
{
    function_params *p = params;

//    realtype E; /* applied potential */
//
//    if (t >= 0 && t <= p->t1)
//        E = p->E0 + p->v * t;
//    if (t > p->t1 && t <= RCONST(2.0) * p->t1)
//        E = p->E0 + RCONST(2.0) * p->v * p->t1 - p->v * t;
//    if (t > RCONST(2.0) * p->t1)
//    {
//        printf("time t is out of range!\n");
//    }

    double result;

    result = p->Ab * gsl_sf_exp( - gsl_pow_2( p->n * GSL_CONST_MKSA_ELECTRON_CHARGE * (E_applied(t, params) - p->E_ref) - p->lamb) / (4.0 * p->lamb * GSL_CONST_MKSA_BOLTZMANN * p->T));
    return result;

}


int main()
{
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
    Ith(y,2) = PR_initial;

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

    /* Call CVodeRootInit to specify the root function g with 2 components */
    flag = CVodeRootInit(cvode_mem, 1, g);

    /* Call CVDense to specify the CVDENSE dense linear solver */
    flag = CVDense(cvode_mem, NEQ);

    /* Set the Jacobian routine to Jac (user-supplied) */
    flag = CVDlsSetDenseJacFn(cvode_mem, Jac);
//    CVodeSetInitStep(cvode_mem, RCONST(1.0e-4));

    /* In loop, call CVode, print results, and test for error.
     Break out of loop when NOUT preset output times have been reached.  */
    printf(" \n Single molecule cyclic voltammetry \n\n");

    iout = 0;  tout = T1;
    while(1) {
        flag = CVode(cvode_mem, tout, y, &t, CV_NORMAL);
        PrintOutput(t, Ith(y,1), Ith(y,2));

        if (flag == CV_ROOT_RETURN) {
            flagr = CVodeGetRootInfo(cvode_mem, rootsfound);
            PrintRootInfo(rootsfound[0]);
        }

        if (check_flag(&flag, "CVode", 1)) break;
        if (flag == CV_SUCCESS) {
            iout++;
            tout += TADD;
//            tout *= TMULT;
        }

        if (iout == NOUT) break;
    }

    /* Print some final statistics */
    PrintFinalStats(cvode_mem);

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

    y1 = Ith(y,1); y2 = Ith(y,2);


    Ith(ydot,1) = - kf(t, &params) * y1 + kb(t, &params) * y2;
    Ith(ydot,2) =   kf(t, &params) * y1 - kb(t, &params) * y2;

    return(0);
}

/*
 * g routine. Compute functions g_i(t,y) for i = 0,1.
 */

static int g(double t, N_Vector y, double *gout, void *user_data)
{
  double y1;

  y1 = Ith(y,1);
  gout [0] = y1 - RCONST(0.5);

  return(0);
}

/*
 * Jacobian routine. Compute J(t,y) = df/dy. *
 */

static int Jac(long int N, double t,
               N_Vector y, N_Vector fy, DlsMat J, void *user_data,
               N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{

    IJth(J,1,1) = - kf(t, &params);
    IJth(J,1,2) =   kb(t, &params);
    IJth(J,2,1) =   kf(t, &params);
    IJth(J,2,2) = - kb(t, &params);

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
  printf("At t = %0.4Le  E = %8.6fV  PO = %14.6Le  PR = %14.6Le\n", t, E_applied(t, &params), y1, y2);
#elif defined(SUNDIALS_DOUBLE_PRECISION)
  printf("At t = %0.4le  E = %8.6fV  PO = %14.6le  PR = %14.6le\n", t, E_applied(t, &params), y1, y2);
#else
  printf("At t = %0.4e  E = %8.6fV  PO = %14.6e  PR = %14.6e\n", t, E_applied(t, &params), y1, y2);
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
