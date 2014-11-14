/* This file is a part of MolStat, which is distributed under the Creative
   Commons Attribution-NonCommercial 4.0 International Public License.
   MolStat (c) 2014, Northwestern University. */
/**
 * \file counterindex.h
 * \brief Provides tools for indexing the bins of a multi-dimensional
 *    histogram.
 *
 * \author Matthew G.\ Reuter
 * \date October 2014
 */

#ifndef __counter_index_h__
#define __counter_index_h__

#include <vector>

namespace molstat {

/**
 * \brief Index class for the bins of a histogram.
 *
 * This class is essentially an iterator to go through the bins of a
 * histogram. Different numbers of bins in different dimensions are allowed,
 * and the molstat::CounterIndex::arrayOffset() function calculates how this
 * bin (via its index) translates to a one-dimensional array.
 */
class CounterIndex
{
private:
	/// The maximum index in each direction.
	const std::vector<std::size_t> max_index;

	/// The index represented by this counter.
	std::vector<std::size_t> index;

public:
	CounterIndex() = delete;

	/**
	 * \brief Constructor that sets the max index.
	 *
	 * A similar constructor that moves the indices (given a rvalue) could also
	 * be implemented, but the length of these arrays is expected to be small.
	 *
	 * \throw std::invalid_argument if one of the max indices is 0.
	 *
	 * \param[in] max_index_ The maximum indices.
	 */
	CounterIndex(const std::vector<std::size_t> &max_index_);

	/**
	 * \brief Increments the counter to the next index.
	 *
	 * \return The incremented counter.
	 */
	CounterIndex operator++();

	/**
	 * \brief Set the index for a specified dimension.
	 *
	 * \throw std::out_of_range if an invalid dimension is specified, if the new
	 *    index is out of range within its dimension, or if the object is at the
	 *    end.
	 *
	 * \param[in] dim The dimension to acces.
	 * \param[in] val The new index within the dimension.
	 */
	void setIndex(const std::size_t dim, const std::size_t val);

	/**
	 * \brief Access the index for a specified dimension.
	 *
	 * \throw std::out_of_range if an invalid dimension is specified or if the
	 *    molstat::CounterIndex is at its end.
	 *
	 * \param[in] dim The dimension to access.
	 * \return The index in dimension `dim`.
	 */
	std::size_t operator[] (const std::size_t dim) const;

	/// Resets a counter to all 0s.
	void reset();

	/// Returns true if this counter is at the end (no more `++`).
	bool at_end() const;

	/**
	 * \brief Gets an 1-D array offset for this index.
	 *
	 * \throw std::out_of_range if this index is at the end (and is not a valid
	 *    index).
	 *
	 * \return The 1-D array offset.
	 */
	std::size_t arrayOffset() const;
};

} // namespace molstat

#endif