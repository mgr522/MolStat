/* This file is a part of MolStat, which is distributed under the Creative
   Commons Attribution-NonCommercial 4.0 International Public License.
   MolStat (c) 2014, Northwestern University. */
/**
 * \file histogram1d_linear.cc
 * \brief Test suite for a 1D histogram class, linear binning.
 *
 * \test Tests the molstat::Histogram class with linear binning
 *    (molstat::BinLinear), 1D histogram.
 *
 * \author Matthew G.\ Reuter
 * \date October 2014
 */

#include <cassert>
#include <cmath>

#include <general/histogram_tools/counterindex.h>
#include <general/histogram_tools/histogram.h>
#include <general/histogram_tools/bin_linear.h>

using namespace std;

/**
 * \brief Main function for testing the Histogram (1D) class with linear
 *    binning.
 *
 * \param[in] argc The number of command-line arguments.
 * \param[in] argv The command-line arguments.
 * \return Exit status: 0 if the code passes the test, non-zero otherwise.
 */
int main(int argc, char **argv)
{
	shared_ptr<molstat::BinStyle> bstyle{ make_shared<molstat::BinLinear>(5) };
	molstat::Histogram hist(1);
	const double thresh = 1.0e-6;

	// artificially populate the histogram
	hist.add_data({0.00}); // 1
	hist.add_data({1.00}); // 5
	hist.add_data({0.12}); // 1
	hist.add_data({0.87}); // 5
	hist.add_data({0.66}); // 4
	hist.add_data({0.50}); // 3
	hist.add_data({0.92}); // 5
	hist.add_data({0.42}); // 3
	hist.add_data({0.21}); // 2
	hist.add_data({0.18}); // 1
	hist.add_data({0.04}); // 1
	hist.add_data({0.99}); // 5
	hist.add_data({0.77}); // 4

	// try to add data with the wrong dimensionality
	try
	{
		hist.add_data({ 0.4, 0.3 });
		assert(false);
	}
	catch(const invalid_argument &e)
	{
		// should be here
	}

	// try to access an iterator; should fail because we haven't yet binned
	try
	{
		hist.begin();
		assert(false);
	}
	catch(const runtime_error &e)
	{
		// should be here
	}

	// make a dummy iterator and verify that we can't access coordinates
	// or bin counts before binning
	try
	{
		molstat::CounterIndex ci{ {2} };
		hist.getCoordinates(ci);
		assert(false);
	}
	catch(const runtime_error &e)
	{
		// should be here
	}

	try
	{
		molstat::CounterIndex ci{ {2} };
		hist.getBinCount(ci);
		assert(false);
	}
	catch(const runtime_error &e)
	{
		// should be here
	}

	// bin the data
	try
	{
		// first, with the wrong number of binning styles
		hist.bin_data({bstyle, bstyle});
		assert(false);
	}
	catch(const invalid_argument &e)
	{
		// should be here
	}

	try
	{
		// second, with a null binning style
		hist.bin_data({ nullptr });
		assert(false);
	}
	catch(const runtime_error &e)
	{
		// should be here
	}

	// now, actually bin it
	hist.bin_data({ bstyle });

	// make sure that the data acquisition functions now throw exceptions
	try
	{
		hist.add_data({0.1});
		assert(false);
	}
	catch(const runtime_error &e)
	{
		// should be here
	}

	try
	{
		hist.bin_data({ bstyle });
		assert(false);
	}
	catch(const runtime_error &e)
	{
		// should be here
	}

	// check the bin contents and the iterator class
	molstat::CounterIndex iter = hist.begin();

	// bin 0 (#1 above)
	// the average coordinate is 0.1 and the bin count should be 4
	assert(abs(hist.getCoordinates(iter)[0] - 0.1) < thresh);
	assert(abs(hist.getBinCount(iter) - 4.*bstyle->dmaskdx(0.1)) < thresh);

	++iter;
	// bin 1 (#2 above)
	// the average coordinate is 0.3 and the bin count should be 1
	assert(abs(hist.getCoordinates(iter)[0] - 0.3) < thresh);
	assert(abs(hist.getBinCount(iter) - 1.*bstyle->dmaskdx(0.3)) < thresh);

	++iter;
	// bin 2 (#3 above)
	// the average coordinate is 0.5 and the bin count should be 2
	assert(abs(hist.getCoordinates(iter)[0] - 0.5) < thresh);
	assert(abs(hist.getBinCount(iter) - 2.*bstyle->dmaskdx(0.5)) < thresh);

	++iter;
	// bin 3 (#4 above)
	// the average coordinate is 0.7 and the bin count should be 2
	assert(abs(hist.getCoordinates(iter)[0] - 0.7) < thresh);
	assert(abs(hist.getBinCount(iter) - 2.*bstyle->dmaskdx(0.7)) < thresh);

	++iter;
	// bin 4 (#5 above)
	// the average coordinate is 0.9 and the bin count should be 4
	assert(abs(hist.getCoordinates(iter)[0] - 0.9) < thresh);
	assert(abs(hist.getBinCount(iter) - 4.*bstyle->dmaskdx(0.9)) < thresh);

	// just for sanity
	++iter;
	assert(iter.at_end());

	// -------------------------------
	// attempt to create some bad histograms to check exception throwing.
	molstat::Histogram hist2(1);
	bstyle = make_shared<molstat::BinLinear>(0);

	// fail due to zero bins in a direction
	try
	{
		hist2.bin_data({ bstyle });
		assert(false);
	}
	catch(const runtime_error &e)
	{
		// should be here
	}

	hist2.add_data({0.5});
	bstyle = make_shared<molstat::BinLinear>(2);

	// fail because the range is null and more than 1 bin is requested
	try
	{
		hist2.bin_data({ bstyle });
		assert(false);
	}
	catch(const size_t &dim)
	{
		// should be here
		assert(dim == 0);
	}

	bstyle = make_shared<molstat::BinLinear>(1);
	hist2.bin_data({ bstyle });

	molstat::CounterIndex iter2 { hist2.begin() };
	assert(abs(hist2.getCoordinates(iter2)[0] - 0.5) < thresh);
	assert(abs(hist2.getBinCount(iter2) - 1.*bstyle->dmaskdx(0.5)) < thresh);

	return 0;
}
