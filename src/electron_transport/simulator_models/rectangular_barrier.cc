/* This file is a part of MolStat, which is distributed under the Creative
   Commons Attribution-NonCommercial 4.0 International Public License.

   (c) 2015 Stony Brook University. */

/**
 * \file rectangular_barrier.cc
 * \brief Channel with transmission given by a rectangular barrier.
 *
 * \author Matthew G.\ Reuter
 * \date January 2015
 */

#include "rectangular_barrier.h"
#include <cmath>
#include <gsl/gsl_integration.h>

namespace molstat {
namespace transport {

const std::size_t RectangularBarrier::Index_EF = TransportJunction::Index_EF;
const std::size_t RectangularBarrier::Index_V = TransportJunction::Index_V;
const std::size_t RectangularBarrier::Index_h = 2;
const std::size_t RectangularBarrier::Index_w = 3;

struct my_f_params { double h; double w;};

std::vector<std::string> RectangularBarrier::get_names() const
{
	std::vector<std::string> ret(2);

	// subtract two because we expect 2 parameters from the TransportJunction
	ret[Index_h - 2] = "height";
	ret[Index_w - 2] = "width";

	return ret;
}


double transmissionf (double E, void * p) 
{
	struct my_f_params * params 
    		= (struct my_f_params *)p;
	double h = (params -> h);
	double w = (params -> w);

	return (4. * E * (h - E)) / ((4. * E * (h - E)) + h*h*(sinh(sqrt(h - E) * 5.12317 * w) * sinh(sqrt(h - E) * 5.12317 * w)));
}


double RectangularBarrier::transmission(const double e, const double h,
	const double w)
{
	// sqrt(2m (eV)) / hbar = 5.12317 / nm
	const double sinhval = std::sinh(5.12317 * sqrt(h - e) * w) * h;
	const double intermed = 4. * e * (h - e);

	return intermed / (intermed + sinhval * sinhval);
}

double RectangularBarrier::ZeroBiasG(const std::valarray<double> &params) const
{
	// unpack the parameters
	const double &ef = params[Index_EF];
	const double &h = params[Index_h];
	const double &w = params[Index_w];
	
	return transmission(ef, h, w);
}

double RectangularBarrier::StaticG(const std::valarray<double> &params) const
{
	// unpack the parameters
	const double &ef = params[Index_EF];
	const double &V = params[Index_V];
	const double &h = params[Index_h];
	const double &w = params[Index_w];
	double result, intmin, intmax, error;
	double res, abserr;
	size_t neval;


	gsl_integration_cquad_workspace *ws = 		gsl_integration_cquad_workspace_alloc (1000);
	 
	gsl_function F;
	struct my_f_params p = {h, w};
  	F.function = transmissionf;
  	F.params = &p;

	intmin = (ef-V)/2;
	intmax = (ef+V)/2;
  	gsl_integration_cquad (&F, intmin, intmax, 1e-9, 1e-9,
                        ws, &result, &abserr , &neval); 
	gsl_integration_cquad_workspace_free (ws);
	// conductance represented as a fraction of G_o
	// the integral of the transmission function returns units of eV
	return 2*2.41798926E14/V*result /7.748091723E-5 ;
}

double RectangularBarrier::DistD(const std::valarray<double> &params) const
{
	// unpack the width
	const double &w = params[Index_w];
	
	return w;
}

} // namespace molstat::transport
} // namespace molstat
