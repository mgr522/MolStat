/**
 * \file tests/histogram1d_linear.cc
 * \brief Test suite for the Histogram1D class, linear binning.
 *
 * \author Matthew G.\ Reuter
 * \date May 2014
 */

#include <cstdio>
#include <assert.h>
#include <cmath>

#include <general/histogram_tools/histogram1d.h>
#include <general/histogram_tools/bin_linear.h>

using namespace std;

/**
 * \brief Main function for testing the Histogram1D class with linear binning.
 *
 * \param[in] argc The number of command-line arguments.
 * \param[in] argv The command-line arguments.
 * \return Exit status: 0 if the code passes the test, non-zero otherwise.
 */
int main(int argc, char **argv) {
	shared_ptr<BinStyle> bstyle(make_shared<BinLinear>());
	Histogram1D hist(5, 0., 1., bstyle);
	const double thresh = 1.0e-6;

	// artificially populate the histogram
	hist.add_data(0.00); // 1
	hist.add_data(1.00); // 5 -- outside since 1. is excluded
	hist.add_data(0.12); // 1
	hist.add_data(0.87); // 5
	hist.add_data(0.66); // 4
	hist.add_data(0.50); // 3
	hist.add_data(0.92); // 5
	hist.add_data(0.42); // 3
	hist.add_data(0.21); // 2
	hist.add_data(0.18); // 1
	hist.add_data(0.04); // 1
	hist.add_data(0.99); // 5
	hist.add_data(0.77); // 4

	// check the bin contents and the iterator class
	Histogram1D::const_iterator iter = hist.begin();

	// bin 0 (#1 above)
	// the average coordinate is 0.1 and the bin count should be 4
	assert(abs(iter.get_variable()[0] - 0.1) < thresh);
	assert(abs(iter.get_bin_count() - 4.*bstyle->dudg(0.1)) < thresh);

	++iter;
	// bin 1 (#2 above)
	// the average coordinate is 0.3 and the bin count should be 1
	assert(abs(iter.get_variable()[0] - 0.3) < thresh);
	assert(abs(iter.get_bin_count() - 1.*bstyle->dudg(0.3)) < thresh);

	iter++;
	// bin 2 (#3 above)
	// the average coordinate is 0.5 and the bin count should be 2
	assert(abs(iter.get_variable()[0] - 0.5) < thresh);
	assert(abs(iter.get_bin_count() - 2.*bstyle->dudg(0.5)) < thresh);

	++iter;
	// bin 3 (#4 above)
	// the average coordinate is 0.7 and the bin count should be 2
	assert(abs(iter.get_variable()[0] - 0.7) < thresh);
	assert(abs(iter.get_bin_count() - 2.*bstyle->dudg(0.7)) < thresh);

	iter++;
	// bin 4 (#5 above)
	// the average coordinate is 0.9 and the bin count should be 3
	assert(abs(iter.get_variable()[0] - 0.9) < thresh);
	assert(abs(iter.get_bin_count() - 3.*bstyle->dudg(0.9)) < thresh);

	// just for sanity
	++iter;
	assert(iter == hist.end());

	return 0;
}
