/* This file is a part of MolStat, which is distributed under the Creative
   Commons Attribution-NonCommercial 4.0 International Public License.
   MolStat (c) 2014, Northwestern University. */
/**
 * \file histogram1d.cc
 * \brief Implementation of an interface for the GSL histogram functions.
 *
 * \author Matthew G.\ Reuter
 * \date May 2014
 */

#include "histogram1d.h"

using namespace std;

Histogram1D::Histogram1D(const std::size_t nbin, const double minval,
	const double maxval, const std::shared_ptr<const BinStyle> bstyle)
	: Histogram<1>(bstyle),
	  hist(gsl_histogram_alloc(nbin), &gsl_histogram_free) {

	gsl_histogram_set_ranges_uniform(hist.get(), bstyle->gmask(minval),
		bstyle->gmask(maxval));
}

void Histogram1D::add_data(const double v) {
	gsl_histogram_increment(hist.get(), bstyle->gmask(v));
}

void Histogram1D::add_data(const std::array<double, 1> &v) {
	add_data(v[0]);
}

Histogram1D::const_iterator Histogram1D::begin() const {
	return const_iterator(hist, bstyle);
}

Histogram1D::const_iterator Histogram1D::end() const {
	const_iterator ret(hist, bstyle);

	// move the bin index to the end
	ret.bin[0] = gsl_histogram_bins(hist.get());

	return ret;
}

void Histogram1D::const_iterator::next_bin() {
	// increment the bin
	++bin[0];
}

void Histogram1D::const_iterator::set_output() {
	// if we're at the end, there's no data to set
	if(bin[0] < gsl_histogram_bins(hist.get())) {
		double upper, lower;

		gsl_histogram_get_range(hist.get(), bin[0], &lower, &upper);
		val[0] = 0.5*(bstyle->invgmask(upper) + bstyle->invgmask(lower));

		bincount = gsl_histogram_get(hist.get(), bin[0])
			* bstyle->dudg(val[0]);
	}
}

Histogram1D::const_iterator::const_iterator(
	const std::shared_ptr<const gsl_histogram> h,
	const std::shared_ptr<const BinStyle> bstyle)
	: Histogram<1>::const_iterator(bstyle), hist(h) {

	set_output();
}

Histogram1D::const_iterator Histogram1D::const_iterator::operator++() {
	const_iterator ret = *this;
	next_bin();
	set_output();
	return ret;
}

Histogram1D::const_iterator Histogram1D::const_iterator::operator++(int) {
	next_bin();
	set_output();
	return *this;
}

bool Histogram1D::const_iterator::operator== (const const_iterator &rhs) const
{
	if(hist == rhs.hist) {
		return (bin == rhs.bin);
	}
		return false;
}

bool Histogram1D::const_iterator::operator!= (const const_iterator &rhs) const
{
	return !operator==(rhs);
}
