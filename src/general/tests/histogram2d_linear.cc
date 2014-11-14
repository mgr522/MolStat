/* This file is a part of MolStat, which is distributed under the Creative
   Commons Attribution-NonCommercial 4.0 International Public License.
   MolStat (c) 2014, Northwestern University. */
/**
 * \internal
 * \file histogram2d_linear.cc
 * \brief Test suite for a 2D histogram class (linear binning on both axes).
 *
 * \test Tests the molstat::Histogram class with linear binning
 *    (molstat::BinLinear) on both axes, 2D histogram.
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

using namespace std;

/**
 * \internal
 * \brief Main function for testing the Histogram (2D) class with linear
 *    binning on both axes.
 *
 * \param[in] argc The number of command-line arguments.
 * \param[in] argv The command-line arguments.
 * \return Exit status: 0 if the code passes the test, non-zero otherwise.
 * \endinternal
 */
int main(int argc, char **argv)
{
	shared_ptr<molstat::BinStyle> bstyle(make_shared<molstat::BinLinear>(2));
	molstat::Histogram hist(2);
	const double thresh = 1.0e-6;

	// artificially populate the histogram
	hist.add_data({{0.4, 0.4}}); // 1
	hist.add_data({{0.3, 0.7}}); // 2
	hist.add_data({{0.4, 0.0}}); // 1
	hist.add_data({{1.0, 0.3}}); // 3
	hist.add_data({{0.1, 0.8}}); // 2
	hist.add_data({{0.6, 0.1}}); // 3
	hist.add_data({{0.2, 0.2}}); // 1
	hist.add_data({{0.3, 0.0}}); // 1
	hist.add_data({{0.7, 1.0}}); // 4
	hist.add_data({{0.0, 0.8}}); // 2

	// bin the data
	hist.bin_data({ bstyle, bstyle });

	// check the bin contents and the iterator class
	molstat::CounterIndex iter = hist.begin();

	// bin 0, 0 (#1 above)
	// the average is (0.25, 0.25) and the bin count should be 4
	assert(abs(hist.getCoordinates(iter)[0] - 0.25) < thresh);
	assert(abs(hist.getCoordinates(iter)[1] - 0.25) < thresh);
	assert(abs(hist.getBinCount(iter) -
		4. * bstyle->dmaskdx(0.25) * bstyle->dmaskdx(0.25)) < thresh);

	++iter;
	// bin 1, 0 (#3 above)
	// the average is (0.75, 0.25) and the bin count should be 2
	assert(abs(hist.getCoordinates(iter)[0] - 0.75) < thresh);
	assert(abs(hist.getCoordinates(iter)[1] - 0.25) < thresh);
	assert(abs(hist.getBinCount(iter) -
		2. * bstyle->dmaskdx(0.75) * bstyle->dmaskdx(0.25)) < thresh);

	++iter;
	// bin 0, 1 (#2 above)
	// the average is (0.25, 0.75) and the bin count should be 3
	assert(abs(hist.getCoordinates(iter)[0] - 0.25) < thresh);
	assert(abs(hist.getCoordinates(iter)[1] - 0.75) < thresh);
	assert(abs(hist.getBinCount(iter) -
		3. * bstyle->dmaskdx(0.25) * bstyle->dmaskdx(0.75)) < thresh);

	++iter;
	// bin 1, 1 (#4 above)
	// the average is (0.75, 0.75) and the bin count should be 1
	assert(abs(hist.getCoordinates(iter)[0] - 0.75) < thresh);
	assert(abs(hist.getCoordinates(iter)[1] - 0.75) < thresh);
	assert(abs(hist.getBinCount(iter) -
		1. * bstyle->dmaskdx(0.75) * bstyle->dmaskdx(0.75)) < thresh);

	// just for sanity
	++iter;
	assert(iter.at_end());

	return 0;
}
