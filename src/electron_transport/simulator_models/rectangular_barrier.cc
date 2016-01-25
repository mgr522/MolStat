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

double RectangularBarrier::DispW(const std::valarray<double> &params) const
{
	// unpack the width
	const double &w = params[Index_w];
	
	return w;
}

#if HAVE_GSL
struct my_f_params { double h; double w; };

/**
 * \brief GSL integrator binding for the tranmission function.
 *
 * \param[in] E The energy (integration variable).
 * \param[in] p The parameters (height and width of the barrier).
 * \return The transmission probability through the barrier.
 */
static double rb_trans_gsl(double E, void *p) 
{
	const struct my_f_params * params = (struct my_f_params*)p;
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
	gsl_integration_cquad_workspace *ws = 
		gsl_integration_cquad_workspace_alloc(1000);
	 
	gsl_function F;
	struct my_f_params p {h, w};
  F.function = &rb_trans_gsl;
  F.params = &p;

  // perform the integration
	intmin = ef - 0.5*V;
	intmax = ef + 0.5*V;
  gsl_integration_cquad(&F, intmin, intmax, 1e-9, 1e-9,
                        ws, &result, &abserr, &neval); 

  // free resources
	gsl_integration_cquad_workspace_free(ws);

	return result / V;
}
#endif

} // namespace molstat::transport
} // namespace molstat
