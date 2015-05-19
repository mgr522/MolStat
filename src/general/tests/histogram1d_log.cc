/* This file is a part of MolStat, which is distributed under the Creative
   Commons Attribution-NonCommercial 4.0 International Public License.

   (c) 2014 Northwestern University. */

/**
 * \file histogram1d_log.cc
 * \brief Test suite for the 1D histogram1D, logarithmic binning.
 *
 * \test Tests the molstat::Histogram class with logarithmic binning
 *    (molstat::BinLog), 1D histogram.
 *
 * \author Matthew G.\ Reuter
 * \date October 2014
 */

#include <cassert>
#include <cmath>

#include <general/histogram_tools/counterindex.h>
#include <general/histogram_tools/histogram.h>
#include <general/histogram_tools/bin_log.h>

using namespace std;

/**
 * \brief Main function for testing the Histogram (1D) class with logarithmic
 *    binning.
 *
 * \param[in] argc The number of command-line arguments.
 * \param[in] argv The command-line arguments.
 * \return Exit status: 0 if the code passes the test, non-zero otherwise.
 */
int main(int argc, char **argv)
{
	shared_ptr<molstat::BinStyle> bstyle(make_shared<molstat::BinLog>(5, 10.));
	molstat::Histogram hist(1);
	const double thresh = 1.0e-6;

	// artificially populate the histogram
	hist.add_data({4.e-3}); // 3
	hist.add_data({5.e-1}); // 5
	hist.add_data({1.e-5}); // 1
	hist.add_data({8.e-5}); // 1
	hist.add_data({6.e-3}); // 3
	hist.add_data({7.e-2}); // 4
	hist.add_data({4.e-3}); // 3
	hist.add_data({2.e-5}); // 1
	hist.add_data({1.e+0}); // 5
	hist.add_data({3.e-1}); // 5
	hist.add_data({4.e-3}); // 3
	hist.add_data({6.e-5}); // 1
	hist.add_data({1.e-2}); // 4

	// bin the data
	hist.bin_data({ bstyle });

	// check the bin contents and the iterator class
	molstat::CounterIndex iter = hist.begin();

	// bin 0 (#1 above)
	// the average coordinate is 5.5e-5 and the bin count should be 4
	assert(abs(hist.getCoordinates(iter)[0] - 5.5e-5) < thresh);
	assert(abs(hist.getBinCount(iter) - 4.*bstyle->dmaskdx(5.5e-5)) < thresh);

	++iter;
	// bin 1 (#2 above)
	// the average coordinate is 5.5e-4 and the bin count should be 0
	assert(abs(hist.getCoordinates(iter)[0] - 5.5e-4) < thresh);
	assert(abs(hist.getBinCount(iter) - 0.*bstyle->dmaskdx(5.5e-4)) < thresh);

	++iter;
	// bin 2 (#3 above)
	// the average coordinate is 5.5e-3 and the bin count should be 4
	assert(abs(hist.getCoordinates(iter)[0] - 5.5e-3) < thresh);
	assert(abs(hist.getBinCount(iter) - 4.*bstyle->dmaskdx(5.5e-3)) < thresh);

	++iter;
	// bin 3 (#4 above)
	// the average coordinate is 5.5e-2 and the bin count should be 2
	assert(abs(hist.getCoordinates(iter)[0] - 5.5e-2) < thresh);
	assert(abs(hist.getBinCount(iter) - 2.*bstyle->dmaskdx(5.5e-2)) < thresh);

	++iter;
	// bin 4 (#5 above)
	// the average coordinate is 5.5e-1 and the bin count should be 3
	assert(abs(hist.getCoordinates(iter)[0] - 5.5e-1) < thresh);
	assert(abs(hist.getBinCount(iter) - 3.*bstyle->dmaskdx(5.5e-1)) < thresh);

	// just for sanity
	++iter;
	assert(iter.at_end());

	return 0;
}
