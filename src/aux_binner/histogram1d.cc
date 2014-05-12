/**
 * \file aux_binner/histogram1d.cc
 * \brief Implementation of an interface for the GSL histogram functions.
 *
 * \author Matthew G.\ Reuter
 * \date May 2014
 */

#include "histogram1d.h"
#include <limits>

using namespace std;

Histogram1D::Histogram1D()
	: minval(numeric_limits<double>::max()),
	  maxval(-numeric_limits<double>::max()), hist(nullptr) {
}

void Histogram1D::add_data(const double v) {
	// do we store the data or add it to the histogram?
	if(hist == nullptr) {
		// no histogram, yet... store the data in the queue
		if(v < minval)
			minval = v;
		if(v > maxval)
			maxval = v;

		data.emplace(v);
	}
	else {
		gsl_histogram_increment(hist.get(), v);
	}
}

void Histogram1D::bin(const size_t n) {
	// allocate the GSL histogram object
	// this will destroy (possibly -- shared_ptr) an existing histogram
	hist.reset(gsl_histogram_alloc(n), &gsl_histogram_free);

	// first set the ranges of the bins
	gsl_histogram_set_ranges_uniform(hist.get(), minval, maxval);

	// now add all of the data
	while(!data.empty()) {
		double thisdata = data.front();
		data.pop();
		gsl_histogram_increment(hist.get(), thisdata);
	}
}

Histogram1D::const_iterator Histogram1D::begin() const {
	return const_iterator(hist);
}

Histogram1D::const_iterator Histogram1D::end() const {
	const_iterator ret(hist);

	if(hist != nullptr) {
		// set the bins to the end of the line
		ret.bin = gsl_histogram_bins(hist.get());
	}

	return ret;
}

void Histogram1D::const_iterator::next_bin() {
	// if there is no histogram, there's nothing to iterate over
	if(hist != nullptr) {
		// increment the bin
		++bin;
	}
}

void Histogram1D::const_iterator::set_output() {
	// if there is no histogram, there's no data to set
	// also, if there is a histogram but we're at the end, there's no data to
	// set
	if(hist != nullptr && bin < gsl_histogram_bins(hist.get())) {
		double upper, lower;

		bincount_ = gsl_histogram_get(hist.get(), bin);

		gsl_histogram_get_range(hist.get(), bin, &lower, &upper);
		val_ = 0.5*(upper + lower);
	}
}

Histogram1D::const_iterator::const_iterator()
	: hist(nullptr), bin(0) {
}

Histogram1D::const_iterator::const_iterator(
	const std::shared_ptr<const gsl_histogram> h)
	: hist(h), bin(0) {

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

double Histogram1D::const_iterator::variable() const {
	return val_;
}

double Histogram1D::const_iterator::bin_count() const {
	return bincount_;
}

bool Histogram1D::const_iterator::operator== (const const_iterator &rhs) const
{
	if(hist == rhs.hist) {
		if(hist == nullptr)
			return true;

		return (bin == rhs.bin);
	}
		return false;
}

bool Histogram1D::const_iterator::operator!= (const const_iterator &rhs) const
{
	return !operator==(rhs);
}
