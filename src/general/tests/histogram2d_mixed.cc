/* This file is a part of MolStat, which is distributed under the Creative
   Commons Attribution-NonCommercial 4.0 International Public License.
   MolStat (c) 2014, Northwestern University. */
/**
 * \internal
 * \file histogram2d_mixed.cc
 * \brief Test suite for a 2D Histogram class (mixed linear and logarithmic
 *    binning).
 *
 * \test Tests the molstat::Histogram class with mixed linear
 *    (molstat::BinLinear) and logarithmic binning (molstat::BinLog).
 *
 * \author Matthew G.\ Reuter
 * \date October 2014
 * \endinternal
 */

#include <cassert>
#include <cmath>

#include <general/histogram_tools/counterindex.h>
#include <general/histogram_tools/histogram.h>
#include <general/histogram_tools/bin_linear.h>
#include <general/histogram_tools/bin_log.h>

using namespace std;

/**
 * \internal
 * \brief Main function for testing the Histogram (2D) class with mixed linear
 *    and logarithmic binning.
 *
 * \param[in] argc The number of command-line arguments.
 * \param[in] argv The command-line arguments.
 * \return Exit status: 0 if the code passes the test, non-zero otherwise.
 * \endinternal
 */
int main(int argc, char **argv) {
	shared_ptr<molstat::BinStyle> blog(make_shared<molstat::BinLog>(4, 10.)),
	                              blin(make_shared<molstat::BinLinear>(2));
	molstat::Histogram hist(2);
	const double thresh = 1.0e-6;

	// artificially populate the histogram
	hist.add_data({{0.4, 4.e-3}}); // 1
	hist.add_data({{0.0, 2.e-1}}); // 3
	hist.add_data({{0.6, 8.e-4}}); // 4
	hist.add_data({{0.2, 3.e-2}}); // 2
	hist.add_data({{0.9, 8.e-3}}); // 5
	hist.add_data({{0.8, 4.e-1}}); // 7
	hist.add_data({{0.1, 1.e-4}}); // 0
	hist.add_data({{0.4, 5.e-2}}); // 2
	hist.add_data({{0.2, 4.e-3}}); // 1
	hist.add_data({{1.0, 7.e-4}}); // 4
	hist.add_data({{0.8, 8.e-2}}); // 6
	hist.add_data({{0.2, 4.e-4}}); // 0
	hist.add_data({{0.3, 4.e-2}}); // 2
	hist.add_data({{0.4, 2.e-1}}); // 3
	hist.add_data({{0.6, 9.e-4}}); // 4
	hist.add_data({{0.7, 2.e-2}}); // 6
	hist.add_data({{0.8, 4.e-3}}); // 5
	hist.add_data({{0.2, 1.e+0}}); // 3
	hist.add_data({{0.1, 6.e-3}}); // 1
	hist.add_data({{0.9, 4.e-2}}); // 6
	hist.add_data({{0.6, 3.e-1}}); // 7
	hist.add_data({{0.2, 4.e-4}}); // 0
	hist.add_data({{0.7, 9.e-1}}); // 7
	hist.add_data({{0.8, 5.e-4}}); // 4
	hist.add_data({{0.4, 4.e-3}}); // 1

	// bin the data
	hist.bin_data({ blin, blog });

	// check the bin contents and the iterator class
	molstat::CounterIndex iter = hist.begin();

	// bin 0, 0 (#0 above)
	// the average is (0.25, 5.5e-4) and the bin count should be 3
	assert(abs(hist.getCoordinates(iter)[0] - 0.25) < thresh);
	assert(abs(hist.getCoordinates(iter)[1] - 5.5e-4) < thresh);
	assert(abs(hist.getBinCount(iter) -
		3. * blin->dmaskdx(0.25) * blog->dmaskdx(5.5e-4)) < thresh);

	++iter;
	// bin 1, 0 (#4 above)
	// the average is (0.75, 5.5e-4) and the bin count should be 4
	assert(abs(hist.getCoordinates(iter)[0] - 0.75) < thresh);
	assert(abs(hist.getCoordinates(iter)[1] - 5.5e-4) < thresh);
	assert(abs(hist.getBinCount(iter) -
		4. * blin->dmaskdx(0.75) * blog->dmaskdx(5.5e-4)) < thresh);

	++iter;
	// bin 0, 1 (#1 above)
	// the average is (0.25, 5.5e-3) and the bin count should be 4
	assert(abs(hist.getCoordinates(iter)[0] - 0.25) < thresh);
	assert(abs(hist.getCoordinates(iter)[1] - 5.5e-3) < thresh);
	assert(abs(hist.getBinCount(iter) -
		4. * blin->dmaskdx(0.25) * blog->dmaskdx(5.5e-3)) < thresh);

	++iter;
	// bin 1, 1 (#5 above)
	// the average is (0.75, 5.5e-3) and the bin count should be 2
	assert(abs(hist.getCoordinates(iter)[0] - 0.75) < thresh);
	assert(abs(hist.getCoordinates(iter)[1] - 5.5e-3) < thresh);
	assert(abs(hist.getBinCount(iter) -
		2. * blin->dmaskdx(0.75) * blog->dmaskdx(5.5e-3)) < thresh);

	++iter;
	// bin 0, 2 (#2 above)
	// the average is (0.25, 5.5e-2) and the bin count should be 3
	assert(abs(hist.getCoordinates(iter)[0] - 0.25) < thresh);
	assert(abs(hist.getCoordinates(iter)[1] - 5.5e-2) < thresh);
	assert(abs(hist.getBinCount(iter) -
		3. * blin->dmaskdx(0.25) * blog->dmaskdx(5.5e-2)) < thresh);

	++iter;
	// bin 1, 2 (#6 above)
	// the average is (0.75, 5.5e-2) and the bin count should be 3
	assert(abs(hist.getCoordinates(iter)[0] - 0.75) < thresh);
	assert(abs(hist.getCoordinates(iter)[1] - 5.5e-2) < thresh);
	assert(abs(hist.getBinCount(iter) -
		3. * blin->dmaskdx(0.75) * blog->dmaskdx(5.5e-2)) < thresh);

	++iter;
	// bin 0, 3 (#3 above)
	// the average is (0.25, 5.5e-1) and the bin count should be 3
	assert(abs(hist.getCoordinates(iter)[0] - 0.25) < thresh);
	assert(abs(hist.getCoordinates(iter)[1] - 5.5e-1) < thresh);
	assert(abs(hist.getBinCount(iter) -
		3. * blin->dmaskdx(0.25) * blog->dmaskdx(5.5e-1)) < thresh);

	++iter;
	// bin 1, 3 (#7 above)
	// the average is (0.75, 5.5e-1) and the bin count should be 3
	assert(abs(hist.getCoordinates(iter)[0] - 0.75) < thresh);
	assert(abs(hist.getCoordinates(iter)[1] - 5.5e-1) < thresh);
	assert(abs(hist.getBinCount(iter) -
		3. * blin->dmaskdx(0.75) * blog->dmaskdx(5.5e-1)) < thresh);

	// just for sanity
	++iter;
	assert(iter.at_end());

	return 0;
}
