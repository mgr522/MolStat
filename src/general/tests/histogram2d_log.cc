/* This file is a part of MolStat, which is distributed under the Creative
   Commons Attribution-NonCommercial 4.0 International Public License.
   MolStat (c) 2014, Northwestern University. */
/**
 * \internal
 * \file histogram2d_log.cc
 * \brief Test suite for a 2D Histogram class (logarithmic binning).
 *
 * \test Tests the molstat::Histogram class with logarithmic binning
 *    (molstat::BinLog).
 *
 * \author Matthew G.\ Reuter
 * \date October 2014
 * \endinternal
 */

#include <cassert>
#include <cmath>

#include <general/histogram_tools/counterindex.h>
#include <general/histogram_tools/histogram.h>
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
	shared_ptr<molstat::BinStyle> blog2(make_shared<molstat::BinLog>(2, 10.)),
	                              blog4(make_shared<molstat::BinLog>(4, 10.));
	molstat::Histogram hist(2);
	const double thresh = 1.0e-6;

	// artificially populate the histogram
	hist.add_data({{9.e-3, 4.e-3}}); // 1
	hist.add_data({{1.e-4, 2.e-1}}); // 3
	hist.add_data({{6.e-1, 8.e-4}}); // 4
	hist.add_data({{2.e-3, 3.e-2}}); // 2
	hist.add_data({{9.e-2, 8.e-3}}); // 5
	hist.add_data({{3.e-2, 4.e-1}}); // 7
	hist.add_data({{2.e-3, 1.e-4}}); // 0
	hist.add_data({{8.e-3, 5.e-2}}); // 2
	hist.add_data({{5.e-4, 4.e-3}}); // 1
	hist.add_data({{1.e+0, 7.e-4}}); // 4
	hist.add_data({{4.e-2, 8.e-2}}); // 6
	hist.add_data({{2.e-3, 4.e-4}}); // 0
	hist.add_data({{8.e-3, 4.e-2}}); // 2
	hist.add_data({{2.e-4, 2.e-1}}); // 3
	hist.add_data({{2.e-2, 9.e-4}}); // 4
	hist.add_data({{9.e-1, 2.e-2}}); // 6
	hist.add_data({{3.e-1, 4.e-3}}); // 5
	hist.add_data({{7.e-3, 1.e+0}}); // 3
	hist.add_data({{6.e-4, 6.e-3}}); // 1
	hist.add_data({{8.e-1, 4.e-2}}); // 6
	hist.add_data({{3.e-2, 3.e-1}}); // 7
	hist.add_data({{4.e-3, 4.e-4}}); // 0
	hist.add_data({{8.e-2, 9.e-1}}); // 7
	hist.add_data({{9.e-1, 5.e-4}}); // 4
	hist.add_data({{8.e-3, 4.e-3}}); // 1

	// bin the data
	hist.bin_data({ blog2, blog4 });

	// check the bin contents and the iterator class
	molstat::CounterIndex iter = hist.begin();

	// bin 0, 0 (#0 above)
	// the average is (5.05e-3, 5.5e-4) and the bin count should be 3
	assert(abs(hist.getCoordinates(iter)[0] - 5.05e-3) < thresh);
	assert(abs(hist.getCoordinates(iter)[1] - 5.5e-4) < thresh);
	assert(abs(hist.getBinCount(iter) -
		3. * blog2->dmaskdx(5.05e-3) * blog4->dmaskdx(5.5e-4)) < thresh);

	++iter;
	// bin 1, 0 (#4 above)
	// the average is (5.05e-1, 5.5e-4) and the bin count should be 4
	assert(abs(hist.getCoordinates(iter)[0] - 5.05e-1) < thresh);
	assert(abs(hist.getCoordinates(iter)[1] - 5.5e-4) < thresh);
	assert(abs(hist.getBinCount(iter) -
		4. * blog2->dmaskdx(5.05e-1) * blog4->dmaskdx(5.5e-4)) < thresh);

	++iter;
	// bin 0, 1 (#1 above)
	// the average is (5.05e-3, 5.5e-3) and the bin count should be 4
	assert(abs(hist.getCoordinates(iter)[0] - 5.05e-3) < thresh);
	assert(abs(hist.getCoordinates(iter)[1] - 5.5e-3) < thresh);
	assert(abs(hist.getBinCount(iter) -
		4. * blog2->dmaskdx(5.05e-3) * blog4->dmaskdx(5.5e-3)) < thresh);

	++iter;
	// bin 1, 1 (#5 above)
	// the average is (5.05e-1, 5.5e-3) and the bin count should be 2
	assert(abs(hist.getCoordinates(iter)[0] - 5.05e-1) < thresh);
	assert(abs(hist.getCoordinates(iter)[1] - 5.5e-3) < thresh);
	assert(abs(hist.getBinCount(iter) -
		2. * blog2->dmaskdx(5.05e-1) * blog4->dmaskdx(5.5e-3)) < thresh);

	++iter;
	// bin 0, 2 (#2 above)
	// the average is (5.05e-3, 5.5e-2) and the bin count should be 3
	assert(abs(hist.getCoordinates(iter)[0] - 5.05e-3) < thresh);
	assert(abs(hist.getCoordinates(iter)[1] - 5.5e-2) < thresh);
	assert(abs(hist.getBinCount(iter) -
		3. * blog2->dmaskdx(5.05e-3) * blog4->dmaskdx(5.5e-2)) < thresh);

	++iter;
	// bin 1, 2 (#6 above)
	// the average is (5.05e-1, 5.5e-2) and the bin count should be 3
	assert(abs(hist.getCoordinates(iter)[0] - 5.05e-1) < thresh);
	assert(abs(hist.getCoordinates(iter)[1] - 5.5e-2) < thresh);
	assert(abs(hist.getBinCount(iter) -
		3. * blog2->dmaskdx(5.05e-1) * blog4->dmaskdx(5.5e-2)) < thresh);

	++iter;
	// bin 0, 3 (#3 above)
	// the average is (5.05e-3, 5.5e-1) and the bin count should be 3
	assert(abs(hist.getCoordinates(iter)[0] - 5.05e-3) < thresh);
	assert(abs(hist.getCoordinates(iter)[1] - 5.5e-1) < thresh);
	assert(abs(hist.getBinCount(iter) -
		3. * blog2->dmaskdx(5.05e-3) * blog4->dmaskdx(5.5e-1)) < thresh);

	++iter;
	// bin 1, 3 (#7 above)
	// the average is (5.05e-1, 5.5e-1) and the bin count should be 3
	assert(abs(hist.getCoordinates(iter)[0] - 5.05e-1) < thresh);
	assert(abs(hist.getCoordinates(iter)[1] - 5.5e-1) < thresh);
	assert(abs(hist.getBinCount(iter) -
		3. * blog2->dmaskdx(5.05e-1) * blog4->dmaskdx(5.5e-1)) < thresh);

	// just for sanity
	++iter;
	assert(iter.at_end());

	return 0;
}
