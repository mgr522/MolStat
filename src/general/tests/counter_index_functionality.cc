/* This file is a part of MolStat, which is distributed under the Creative
   Commons Attribution-NonCommercial 4.0 International Public License.
   MolStat (c) 2014, Northwestern University. */
/**
 * \file counter_index_functionality.cc
 * \brief Test suite for the histogram counter/index class.
 *
 * \test Tests the molstat::CounterIndex class.
 *
 * \author Matthew G.\ Reuter
 * \date October 2014
 */

#include <cassert>
#include <vector>
#include <general/histogram_tools/counterindex.h>

using namespace std;

/**
 * \brief Helper function for checking the increment properties of
 *    molstat::CounterIndex.
 *
 * \throw std::exception if the test fails.
 *
 * \param[in] ci The molstat::CounterIndex object.
 * \param[in] index The indices to check.
 * \param[in] offset The offset for this index.
 */
static void testIndex(const molstat::CounterIndex &ci,
	const vector<size_t> &index, const size_t offset);

/**
 * \brief Main function for testing the molstat::CounterIndex class.
 *
 * \param[in] argc The number of command-line arguments.
 * \param[in] argv The command-line arguments.
 * \return Exit status: 0 if the code passes the test, non-zero otherwise.
 */
int main(int argc, char **argv)
{
	vector<size_t> max{0, 1, 2, 4};

	// try to make a ConstantIndex -- should fail due to the 0 index
	try
	{
		molstat::CounterIndex attempt{ max };
		assert(false);
	}
	catch(const invalid_argument &e)
	{
		// should be here
	}

	// make a real CounterIndex now
	max[0] = 3;
	molstat::CounterIndex ci { max };
	assert(!ci.at_end());

	// check a bad dimension
	try
	{
		ci[5];
		assert(false);
	}
	catch(const std::out_of_range &e)
	{
		// should be here
	}

	try
	{
		ci.setIndex(5, 0); // another bad dimension
		assert(false);
	}
	catch(const std::out_of_range &e)
	{
		// should be here
	}

	// go through the increments and check the behavior
	size_t offset{ 0 };
	testIndex(ci, {0, 0, 0, 0}, offset);
	++ci; ++offset;
	testIndex(ci, {1, 0, 0, 0}, offset);
	++ci; ++offset;
	testIndex(ci, {2, 0, 0, 0}, offset);
	++ci; ++offset;
	testIndex(ci, {0, 0, 1, 0}, offset);
	++ci; ++offset;
	testIndex(ci, {1, 0, 1, 0}, offset);
	++ci; ++offset;
	testIndex(ci, {2, 0, 1, 0}, offset);
	++ci; ++offset;
	testIndex(ci, {0, 0, 0, 1}, offset);
	++ci; ++offset;
	testIndex(ci, {1, 0, 0, 1}, offset);
	++ci; ++offset;
	testIndex(ci, {2, 0, 0, 1}, offset);
	++ci; ++offset;
	testIndex(ci, {0, 0, 1, 1}, offset);
	++ci; ++offset;
	testIndex(ci, {1, 0, 1, 1}, offset);
	++ci; ++offset;
	testIndex(ci, {2, 0, 1, 1}, offset);
	++ci; ++offset;
	testIndex(ci, {0, 0, 0, 2}, offset);
	++ci; ++offset;
	testIndex(ci, {1, 0, 0, 2}, offset);
	++ci; ++offset;
	testIndex(ci, {2, 0, 0, 2}, offset);
	++ci; ++offset;
	testIndex(ci, {0, 0, 1, 2}, offset);
	++ci; ++offset;
	testIndex(ci, {1, 0, 1, 2}, offset);
	++ci; ++offset;
	testIndex(ci, {2, 0, 1, 2}, offset);
	++ci; ++offset;
	testIndex(ci, {0, 0, 0, 3}, offset);
	++ci; ++offset;
	testIndex(ci, {1, 0, 0, 3}, offset);
	++ci; ++offset;
	testIndex(ci, {2, 0, 0, 3}, offset);
	++ci; ++offset;
	testIndex(ci, {0, 0, 1, 3}, offset);
	++ci; ++offset;
	testIndex(ci, {1, 0, 1, 3}, offset);
	++ci; ++offset;
	testIndex(ci, {2, 0, 1, 3}, offset);
	++ci; ++offset;

	assert(ci.at_end());
	++ci;
	assert(ci.at_end());

	// try to do some things that we shouldn't do because we're at the end
	try
	{
		ci.arrayOffset(); // no offset because we're at the end
		assert(false);
	}
	catch(const out_of_range &e)
	{
		// should be here
	}

	try
	{
		ci.setIndex(0, 0); // cannot set when ci is at the end
		assert(false);
	}
	catch(const out_of_range &e)
	{
		// should be here
	}

	try
	{
		ci[0]; // cannot use [] at the end
		assert(false);
	}
	catch(const out_of_range &e)
	{
		// should be here
	}

	// reset to make operations viable again
	ci.reset();
	testIndex(ci, {0, 0, 0, 0}, 0);

	// try to set with an out-of-range index
	try
	{
		ci.setIndex(2, 2);
		assert(false);
	}
	catch(const std::out_of_range &e)
	{
		// should be here
	}

	ci.setIndex(2, 1);
	testIndex(ci, {0, 0, 1, 0}, 3);

	max[0] = 2; // make sure the index doesn't rely on the original size array
	testIndex(ci, {0, 0, 1, 0}, 3);

	return 0;
}

static void testIndex(const molstat::CounterIndex &ci,
	const vector<size_t> &index, const size_t offset)
{
	assert(index.size() == 4);

	for(size_t j = 0; j < index.size(); ++j)
	{
		assert(ci[j] == index[j]);
	}

	assert(ci.arrayOffset() == offset);
}
