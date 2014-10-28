/* This file is a part of MolStat, which is distributed under the Creative
   Commons Attribution-NonCommercial 4.0 International Public License.
   MolStat (c) 2014, Northwestern University. */
/**
 * \file counterindex.cc
 * \brief Implements tools for indexing the bins of a multi-dimensional
 *    histogram.
 *
 * \author Matthew G.\ Reuter
 * \date October 2014
 */

#include "counterindex.h"
#include <string>

namespace molstat {

CounterIndex::CounterIndex(const std::vector<std::size_t> &max_index_)
	: max_index(max_index_), index(max_index_.size() + 1, 0)
	// note: add 1 dimension to index to store the "at end" bit
{
	// check for 0 indices
	for(std::size_t j = 0; j < max_index.size(); ++j)
		if(max_index[j] == 0)
			throw std::invalid_argument("0 max index detected.");
}

CounterIndex CounterIndex::operator++()
{
	if(!at_end()) // nothing to do if we're at the end
	{
		// augment the first index
		++index[0];

		// look for overflows, and keep propagating up the indices if detected
		for(std::size_t j = 0; j < index.size(); ++j)
		{
			if(index[j] == max_index[j])
			{
				index[j] = 0;
				++index[j + 1];
			}
		}
	}

	return *this;
}

void CounterIndex::setIndex(const std::size_t dim, const std::size_t val)
{
	if(at_end())
		throw std::out_of_range("CounterIndex is at the end: cannot set.");

	// make sure the index is in range
	try
	{
		// use at the check the dimension number
		if(val >= max_index.at(dim))
			throw std::out_of_range("Invalid index for dimension " +
				std::to_string(dim) + '.');
	}
	catch(const std::out_of_range &e)
	{
		throw std::out_of_range("Invalid dimension in CounterIndex::setIndex.");
	}

	// don't need to check bounds for the actual set
	index[dim] = val;
}

std::size_t CounterIndex::operator[] (const std::size_t dim) const
{
	if(at_end())
		throw std::out_of_range("CounterIndex is at the end: [] operator " \
			"invalid.");
	
	// at will throw the out_of_range, if necessary
	return index.at(dim);
}

void CounterIndex::reset()
{
	// set all of the index counts to 0
	for(std::size_t j = 0; j < index.size(); ++j)
		index[j] = 0;
}

bool CounterIndex::at_end() const
{
	// the last element of index is a "end" marker.
	return index[index.size() - 1] != 0;
}

std::size_t CounterIndex::arrayOffset() const
{
	if(at_end())
		throw std::out_of_range("No offset for an index at the end.");

	std::size_t ret{ 0 };

	// -2 in the initialization -- remember, there's an extra value for the "at
	// end" index
	for(int j = index.size() - 2; j >= 0; --j)
	{
		ret = max_index[j] * ret + index[j];
	}

	return ret;
}

} // namespace molstat