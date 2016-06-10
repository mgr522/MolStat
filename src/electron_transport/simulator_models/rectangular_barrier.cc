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

#if HAVE_GSL
#include <memory>
#include <gsl/gsl_integration.h>
#endif

namespace molstat {
namespace transport {

const std::size_t RectangularBarrier::Index_EF = TransportJunction::Index_EF;
const std::size_t RectangularBarrier::Index_V = TransportJunction::Index_V;
const std::size_t RectangularBarrier::Index_h = 2;
const std::size_t RectangularBarrier::Index_w = 3;

std::vector<std::string> RectangularBarrier::get_names() const
{
	std::vector<std::string> ret(2);

	// subtract two because we expect 2 parameters from the TransportJunction
	ret[Index_h - 2] = "height";
	ret[Index_w - 2] = "width";

	return ret;
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

double RectangularBarrier::ZeroBiasS(const std::valarray<double> &params) const
{
	// unpack the parameter
	const double &ef = params[Index_EF];
	
	return -1. / ef;
}

double RectangularBarrier::DispW(const std::valarray<double> &params) const
{
	// unpack the width
	const double &w = params[Index_w];
	
	return w;
}

#if HAVE_GSL
double RectangularBarrier::gsl_StaticG_integrand(double E, void *p) 
{
	const StaticG_data *params = (const StaticG_data*)p;
	const double h = params->h;
	const double w = params->w;

	return RectangularBarrier::transmission(E, h, w);
}

double RectangularBarrier::StaticG(const std::valarray<double> &params) const
{
	// unpack the parameters
	const double &ef = params[Index_EF];
	const double &V = params[Index_V];
	const double &h = params[Index_h];
	const double &w = params[Index_w];
	double result, intmin, intmax;
	double abserr;
	size_t neval;

	// set up the GSL workspace
	std::unique_ptr<gsl_integration_cquad_workspace,
			decltype(&gsl_integration_cquad_workspace_free)>
		ws { gsl_integration_cquad_workspace_alloc(1000),
		     &gsl_integration_cquad_workspace_free };
	 
	gsl_function F;
	struct StaticG_data p {h, w};
  F.function = &gsl_StaticG_integrand;
  F.params = &p;

  // perform the integration
	intmin = 0;
	intmax = V;
  gsl_integration_cquad(&F, intmin, intmax, 1e-9, 1e-9,
                        ws.get(), &result, &abserr, &neval); 

	return result / V;
}
#endif

} // namespace molstat::transport
} // namespace molstat
