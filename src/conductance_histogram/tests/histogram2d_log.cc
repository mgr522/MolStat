/**
 * \file tests/histogram2d_log.cc
 * \brief Test suite for the Histogram2D class (logarithmic binning)
 *
 * \author Matthew G.\ Reuter
 * \date May 2014
 */

#include <cstdio>
#include <assert.h>
#include <cmath>

#include <general/histogram_tools/histogram2d.h>
#include <general/histogram_tools/bin_log.h>

using namespace std;

/**
 * \brief Main function for testing the Histogram2D class (logarithmic binning).
 *
 * \param[in] argc The number of command-line arguments.
 * \param[in] argv The command-line arguments.
 * \return Exit status: 0 if the code passes the test, non-zero otherwise.
 */
int main(int argc, char **argv) {
	shared_ptr<BinStyle> bstyle(make_shared<BinLog>(10.));
	Histogram2D hist({{2,4}}, {{0.,1.e-4}}, {{1.,1.}}, bstyle);
	const double thresh = 1.0e-6;

	// artificially populate the histogram
	hist.add_data({{0.4, 4.e-3}}); // 1
	hist.add_data({{0.1, 2.e-1}}); // 3
	hist.add_data({{0.6, 8.e-4}}); // 4
	hist.add_data({{0.2, 3.e-2}}); // 2
	hist.add_data({{0.9, 8.e-3}}); // 5
	hist.add_data({{0.8, 4.e-1}}); // 7
	hist.add_data({{0.1, 2.e-4}}); // 0
	hist.add_data({{0.4, 5.e-2}}); // 2
	hist.add_data({{0.2, 4.e-3}}); // 1
	hist.add_data({{0.9, 7.e-4}}); // 4
	hist.add_data({{0.8, 8.e-2}}); // 6
	hist.add_data({{0.2, 4.e-4}}); // 0
	hist.add_data({{0.3, 4.e-2}}); // 2
	hist.add_data({{0.4, 2.e-1}}); // 3
	hist.add_data({{0.6, 9.e-4}}); // 4
	hist.add_data({{0.7, 2.e-2}}); // 6
	hist.add_data({{0.8, 4.e-3}}); // 5
	hist.add_data({{0.2, 8.e-1}}); // 3
	hist.add_data({{0.1, 6.e-3}}); // 1
	hist.add_data({{0.9, 4.e-2}}); // 6
	hist.add_data({{0.6, 3.e-1}}); // 7
	hist.add_data({{0.2, 4.e-4}}); // 0
	hist.add_data({{0.7, 9.e-1}}); // 7
	hist.add_data({{0.8, 5.e-4}}); // 4
	hist.add_data({{0.4, 4.e-3}}); // 1

	// check the bin contents and the iterator class
	Histogram2D::const_iterator iter = hist.begin();

	// bin 0, 0 (#0 above)
	// the average is (0.25, 5.5e-4) and the bin count should be 3
	assert(abs(iter.get_variable()[0] - 0.25) < thresh);
	assert(abs(iter.get_variable()[1] - 5.5e-4) < thresh);
	assert(abs(iter.get_bin_count() - 3.*bstyle->dudg(5.5e-4)) < thresh);

	++iter;
	// bin 0, 1 (#1 above)
	// the average is (0.25, 5.5e-3) and the bin count should be 4
	assert(abs(iter.get_variable()[0] - 0.25) < thresh);
	assert(abs(iter.get_variable()[1] - 5.5e-3) < thresh);
	assert(abs(iter.get_bin_count() - 4.*bstyle->dudg(5.5e-3)) < thresh);

	iter++;
	// bin 0, 2 (#2 above)
	// the average is (0.25, 5.5e-2) and the bin count should be 3
	assert(abs(iter.get_variable()[0] - 0.25) < thresh);
	assert(abs(iter.get_variable()[1] - 5.5e-2) < thresh);
	assert(abs(iter.get_bin_count() - 3.*bstyle->dudg(5.5e-2)) < thresh);

	++iter;
	// bin 0, 3 (#3 above)
	// the average is (0.25, 5.5e-1) and the bin count should be 3
	assert(abs(iter.get_variable()[0] - 0.25) < thresh);
	assert(abs(iter.get_variable()[1] - 5.5e-1) < thresh);
	assert(abs(iter.get_bin_count() - 3.*bstyle->dudg(5.5e-1)) < thresh);

	iter++;
	// bin 1, 0 (#4 above)
	// the average is (0.75, 5.5e-4) and the bin count should be 4
	assert(abs(iter.get_variable()[0] - 0.75) < thresh);
	assert(abs(iter.get_variable()[1] - 5.5e-4) < thresh);
	assert(abs(iter.get_bin_count() - 4.*bstyle->dudg(5.5e-4)) < thresh);

	++iter;
	// bin 1, 1 (#5 above)
	// the average is (0.75, 5.5e-3) and the bin count should be 2
	assert(abs(iter.get_variable()[0] - 0.75) < thresh);
	assert(abs(iter.get_variable()[1] - 5.5e-3) < thresh);
	assert(abs(iter.get_bin_count() - 2.*bstyle->dudg(5.5e-3)) < thresh);

	iter++;
	// bin 1, 2 (#6 above)
	// the average is (0.75, 5.5e-2) and the bin count should be 3
	assert(abs(iter.get_variable()[0] - 0.75) < thresh);
	assert(abs(iter.get_variable()[1] - 5.5e-2) < thresh);
	assert(abs(iter.get_bin_count() - 3.*bstyle->dudg(5.5e-2)) < thresh);

	++iter;
	// bin 1, 3 (#7 above)
	// the average is (0.75, 5.5e-1) and the bin count should be 3
	assert(abs(iter.get_variable()[0] - 0.75) < thresh);
	assert(abs(iter.get_variable()[1] - 5.5e-1) < thresh);
	assert(abs(iter.get_bin_count() - 3.*bstyle->dudg(5.5e-1)) < thresh);

	// just for sanity
	++iter;
	assert(iter == hist.end());

	return 0;
}
