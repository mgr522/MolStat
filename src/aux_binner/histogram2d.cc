/**
 * \file histogram2d.cc
 * \brief Implementation of a class for making 2D histograms from data.
 *
 * \author Matthew G.\ Reuter
 * \date May 2014
 */

#include "histogram2d.h"
#include <limits>

using namespace std;

Histogram2D::Histogram2D()
	: mins(numeric_limits<double>::max(), numeric_limits<double>::max()),
	  maxs(-numeric_limits<double>::max(), -numeric_limits<double>::max()),
	  hist(nullptr) {
}

void Histogram2D::add_data(const double v1, const double v2) {
	// do we store the data or add it to the histogram?
	if(hist == nullptr) {
		// no histogram, yet... store the data in the queue
		if(v1 < mins.first)
			mins.first = v1;
		if(v1 > maxs.first)
			maxs.first = v1;
		if(v2 < mins.second)
			mins.second = v2;
		if(v2 > maxs.second)
			maxs.second = v2;

		data.emplace(v1, v2);
	}
	else {
		gsl_histogram2d_increment(hist.get(), v1, v2);
	}
}

void Histogram2D::bin(const size_t n1, const size_t n2) {
	// allocate the GSL histogram object
	// this will destroy (possibly -- shared_ptr) an existing histogram
	hist.reset(gsl_histogram2d_alloc(n1, n2), &gsl_histogram2d_free);

	// first set the ranges of the bins
	gsl_histogram2d_set_ranges_uniform(hist.get(), mins.first, maxs.first,
		mins.second, maxs.second);

	// now add all of the data
	while(!data.empty()) {
		pair<double, double> thisdata = data.front();
		data.pop();
		gsl_histogram2d_increment(hist.get(), thisdata.first, thisdata.second);
	}
}

Histogram2D::const_iterator Histogram2D::begin() const {
	return const_iterator(hist);
}

Histogram2D::const_iterator Histogram2D::end() const {
	const_iterator ret(hist);

	if(hist != nullptr) {
		// set the bins to the end of the line
		ret.bin1 = gsl_histogram2d_nx(hist.get());
		ret.bin2 = gsl_histogram2d_ny(hist.get());
	}

	return ret;
}

void Histogram2D::const_iterator::next_bin() {
	// if there is no histogram, there's nothing to iterate over
	if(hist != nullptr) {
		// increment the bin for the second variable
		++bin2;

		// check if we need to "carry" into the bin for the first variable
		if(bin2 >= gsl_histogram2d_ny(hist.get())) {
			++bin1;

			// check if this is the end; if it ISN'T, reset bin2
			if(bin1 < gsl_histogram2d_nx(hist.get()))
				bin2 = 0;
			// now, if this is the end, bin1 and bin2 equal the max values
		}
	}
}

void Histogram2D::const_iterator::set_output() {
	// if there is no histogram, there's no data to set
	// also, if there is a histogram but we're at the end, there's no data to
	// set
	if(hist != nullptr && bin1 < gsl_histogram2d_nx(hist.get()) &&
		bin2 < gsl_histogram2d_ny(hist.get())) {

		double upper, lower;

		bincount_ = gsl_histogram2d_get(hist.get(), bin1, bin2);

		gsl_histogram2d_get_xrange(hist.get(), bin1, &lower, &upper);
		val1_ = 0.5*(upper + lower);

		gsl_histogram2d_get_yrange(hist.get(), bin2, &lower, &upper);
		val2_ = 0.5*(upper + lower);
	}
}

Histogram2D::const_iterator::const_iterator()
	: hist(nullptr), bin1(0), bin2(0) {
}

Histogram2D::const_iterator::const_iterator(
	const std::shared_ptr<const gsl_histogram2d> h)
	: hist(h), bin1(0), bin2(0) {

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

double Histogram2D::const_iterator::variable1() const {
	return val1_;
}

double Histogram2D::const_iterator::variable2() const {
	return val2_;
}

double Histogram2D::const_iterator::bin_count() const {
	return bincount_;
}

bool Histogram2D::const_iterator::operator== (const const_iterator &rhs) const
{
	if(hist == rhs.hist) {
		if(hist == nullptr)
			return true;

		return (bin1 == rhs.bin1 && bin2 == rhs.bin2);
	}
		return false;
}

bool Histogram2D::const_iterator::operator!= (const const_iterator &rhs) const
{
	return !operator==(rhs);
}
