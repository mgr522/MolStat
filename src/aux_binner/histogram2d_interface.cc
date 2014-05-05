/**
 * \file histogram2d_interface.cc
 * \brief Implementation of a class for making 2D histograms from data.
 *
 * \author Matthew G.\ Reuter
 * \date May 2014
 */

#include "histogram2d_interface.h"
#include <limits>

using namespace std;

Histogram2D::Histogram2D(const int n1, const int n2)
	: mins(numeric_limits<double>::max(), numeric_limits<double>::max()),
	  maxs(-numeric_limits<double>::max(), -numeric_limits<double>::max()),
	  hist(gsl_histogram2d_alloc(n1, n2), &gsl_histogram2d_free) {
}

void Histogram2D::add_data(const double v1, const double v2) {
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

void Histogram2D::bin() {
}

Histogram2D::const_iterator Histogram2D::begin() const {
	return const_iterator(hist);
}

Histogram2D::const_iterator Histogram2D::end() const {
	return const_iterator(hist);
}

void Histogram2D::const_iterator::next_bin() {
}

void Histogram2D::const_iterator::set_output() {
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
	return true;
}

bool Histogram2D::const_iterator::operator!= (const const_iterator &rhs) const
{
	return !operator==(rhs);
}
