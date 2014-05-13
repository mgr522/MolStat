/**
 * \file tests/histogram1d_log.cc
 * \brief Test suite for the Histogram1D class, logarithmic binning.
 *
 * \author Matthew G.\ Reuter
 * \date May 2014
 */

#include <cstdio>
#include <assert.h>
#include <cmath>
#include "../aux_binner/histogram1d.h"
#include "../aux_binner/bin_log.h"

using namespace std;

/**
 * \brief Main function for testing the Histogram1D class with logarithmic
 *    binning.
 *
 * \param[in] argc The number of command-line arguments.
 * \param[in] argv The command-line arguments.
 * \return Exit status: 0 if the code passes the test, non-zero otherwise.
 */
int main(int argc, char **argv) {
	shared_ptr<BinStyle> bstyle(make_shared<BinLog>(10.));
	Histogram1D hist(5, 1.e-5, 1., bstyle);
	const double thresh = 1.0e-6;

	// artificially populate the histogram
	hist.add_data(4.e-3); // 3
	hist.add_data(5.e-1); // 5
	hist.add_data(2.e-5); // 1
	hist.add_data(8.e-5); // 1
	hist.add_data(6.e-3); // 3
	hist.add_data(7.e-2); // 4
	hist.add_data(4.e-3); // 3
	hist.add_data(2.e-5); // 1
	hist.add_data(9.e-1); // 5
	hist.add_data(3.e-1); // 5
	hist.add_data(4.e-3); // 3
	hist.add_data(6.e-5); // 1
	hist.add_data(1.e-2); // 4

	// check the bin contents and the iterator class
	Histogram1D::const_iterator iter = hist.begin();

	// bin 0 (#1 above)
	// the average coordinate is 5.5e-5 and the bin count should be 4
	assert(abs(iter.get_variable()[0] - 5.5e-5) < thresh);
	assert(abs(iter.get_bin_count() - 4.*bstyle->dudg(5.5e-5)) < thresh);

	++iter;
	// bin 1 (#2 above)
	// the average coordinate is 5.5e-4 and the bin count should be 0
	assert(abs(iter.get_variable()[0] - 5.5e-4) < thresh);
	assert(abs(iter.get_bin_count() - 0.*bstyle->dudg(5.5e-4)) < thresh);

	iter++;
	// bin 2 (#3 above)
	// the average coordinate is 5.5e-3 and the bin count should be 4
	assert(abs(iter.get_variable()[0] - 5.5e-3) < thresh);
	assert(abs(iter.get_bin_count() - 4.*bstyle->dudg(5.5e-3)) < thresh);

	++iter;
	// bin 3 (#4 above)
	// the average coordinate is 5.5e-2 and the bin count should be 2
	assert(abs(iter.get_variable()[0] - 5.5e-2) < thresh);
	assert(abs(iter.get_bin_count() - 2.*bstyle->dudg(5.5e-2)) < thresh);

	iter++;
	// bin 4 (#5 above)
	// the average coordinate is 5.5e-1 and the bin count should be 3
	assert(abs(iter.get_variable()[0] - 5.5e-1) < thresh);
	assert(abs(iter.get_bin_count() - 3.*bstyle->dudg(5.5e-1)) < thresh);

	// just for sanity
	++iter;
	assert(iter == hist.end());

	return 0;
}
