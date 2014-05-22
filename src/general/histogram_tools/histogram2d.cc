/* This file is a part of MolStat, which is distributed under the Creative
   Commons Attribution-NonCommercial 4.0 International Public License.
   MolStat (c) 2014, Northwestern University. */
/**
 * \file histogram2d.cc
 * \brief Implementation of a class for making 2D histograms from data.
 *
 * \author Matthew G.\ Reuter
 * \date May 2014
 */

#include "histogram2d.h"

using namespace std;

Histogram2D::Histogram2D(const std::array<std::size_t, 2> &nbin,
	const std::array<double, 2> &mins, const std::array<double, 2> &maxs,
	const std::shared_ptr<const BinStyle> bstyle)
	: Histogram<2>(bstyle),
	  hist(gsl_histogram2d_alloc(nbin[0], nbin[1]), &gsl_histogram2d_free) {

	gsl_histogram2d_set_ranges_uniform(hist.get(), mins[0], maxs[0],
		bstyle->gmask(mins[1]), bstyle->gmask(maxs[1]));
}

void Histogram2D::add_data(const std::array<double, 2> &v) {
	gsl_histogram2d_increment(hist.get(), v[0], bstyle->gmask(v[1]));
}

Histogram2D::const_iterator Histogram2D::begin() const {
	return const_iterator(hist, bstyle);
}

Histogram2D::const_iterator Histogram2D::end() const {
	const_iterator ret(hist, bstyle);

	// move the bin indices to the end
	ret.bin[0] = gsl_histogram2d_nx(hist.get());
	ret.bin[1] = gsl_histogram2d_ny(hist.get());

	return ret;
}

void Histogram2D::const_iterator::next_bin() {
	// increment the bin for the second variable
	++bin[1];

	// check if we need to "carry" into the bin for the first variable
	if(bin[1] >= gsl_histogram2d_ny(hist.get())) {
		++bin[0];

		// check if this is the end; if it ISN'T, reset bin[1]
		if(bin[0] < gsl_histogram2d_nx(hist.get()))
			bin[1] = 0;
		// now, if this is the end, bin[0] and bin[1] equal the max values
	}
}

void Histogram2D::const_iterator::set_output() {
	// if we're at the end, there's no data to set
	if(bin[0] < gsl_histogram2d_nx(hist.get()) &&
		bin[1] < gsl_histogram2d_ny(hist.get())) {

		double upper, lower;

		gsl_histogram2d_get_xrange(hist.get(), bin[0], &lower, &upper);
		val[0] = 0.5*(upper + lower);

		gsl_histogram2d_get_yrange(hist.get(), bin[1], &lower, &upper);
		val[1] = 0.5*(bstyle->invgmask(upper) + bstyle->invgmask(lower));

		bincount = gsl_histogram2d_get(hist.get(), bin[0], bin[1])
			* bstyle->dudg(val[1]);
	}
}

Histogram2D::const_iterator::const_iterator(
	const std::shared_ptr<const gsl_histogram2d> h,
	const std::shared_ptr<const BinStyle> bstyle)
	: Histogram<2>::const_iterator(bstyle), hist(h) {

	set_output();
}

Histogram2D::const_iterator Histogram2D::const_iterator::operator++() {
	const_iterator ret = *this;
	next_bin();
	set_output();
	return ret;
}

Histogram2D::const_iterator Histogram2D::const_iterator::operator++(int) {
	next_bin();
	set_output();
	return *this;
}

bool Histogram2D::const_iterator::operator== (const const_iterator &rhs) const
{
	if(hist == rhs.hist) {
		return (bin[0] == rhs.bin[0] && bin[1] == rhs.bin[1]);
	}
		return false;
}

bool Histogram2D::const_iterator::operator!= (const const_iterator &rhs) const
{
	return !operator==(rhs);
}
