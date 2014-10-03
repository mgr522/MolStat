#ifndef __EXAMPLES_LIBGSL_PO_H__
#define __EXAMPLES_LIBGSL_PO_H__

#include <sundials/sundials_types.h> /* definition of type realtype */
#include <vector>

static void unpack_parameters(const std::vector<double> &vec, double &e0,
        double &eref, double &lambda, double &af, double &ab, double &v,
        double &n, double &poinitial, double &temperature, double &tlimit);
/* Functions Called by the Solver */

static int f(double t, N_Vector y, N_Vector ydot, void *user_data);

static int g(double t, N_Vector y, double *gout, void *user_data);

static int Jac(long int N, double t,
               N_Vector y, N_Vector fy, DlsMat J, void *user_data,
               N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);

/* Private functions to output results */

static void PrintOutput(double t, double y1, double y2);
static void PrintRootInfo(int root_f1);

/* Private function to print final statistics */

static void PrintFinalStats(void *cvode_mem);

/* Private function to check function return values */

static int check_flag(void *flagvalue, char const *funcname, int opt);

/* the function of forward rate constant */
static double kf(double t, const std::vector<double> &vec);

/* the function of backward rate constant */
static double kb(double t, const std::vector<double> &vec);

/* applied potential */
static double E_applied(double t, const std::vector<double> &vec);

#endif /* !defined __EXAMPLES_LIBGSL_PO_H__ */
