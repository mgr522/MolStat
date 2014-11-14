/* This file is a part of MolStat, which is distributed under the Creative
   Commons Attribution-NonCommercial 4.0 International Public License.
   MolStat (c) 2014, Northwestern University. */
/**
 * \file histogram.h
 * \brief Provides tools for constructing histograms.
 *
 * \author Matthew G.\ Reuter
 * \date October 2014
 */

#ifndef __histogram_h__
#define __histogram_h__

#include <memory>
#include <valarray>
#include <vector>
#include <forward_list>
#include <array>
#include "counterindex.h"

namespace molstat {

// forward declaration
class BinStyle;

/**
 * \brief Class that accumulates data and then bins it into a histogram.
 *
 * This class exists in two states:
 * -# Accumulate and store data that will eventually be binned into a
 *    histogram.
 * -# Bin the data into a histogram using the number of bins and binning style
 *    specified by a molstat::BinStyle object.
 *
 * The reason for this split is that the binning operation needs to know the
 * bounds of the data. As data is added, the class tracks the minimum and
 * maximum values.
 *
 * Histograms of any dimensionality can be constructed.
 */
class Histogram
{
private:
	/// State of the histogram: False if binning has not occurred; true if it
	/// has.
	bool haveBinned;

	/// The dimensionality of the data.
	const std::size_t ndim;

	/// The accumulated data.
	std::forward_list<std::valarray<double>> data;

	/// The minimum and maximum values in each dimension.
	std::vector<std::array<double, 2>> extremes;

	/// The number of bins in each dimension.
	std::vector<std::size_t> nbin_dim;

	/**
	 * \brief The middle of each bin.
	 *
	 * The first index is the dimension, the second index is the bin number.
	 */
	std::vector<std::vector<double>> bin_value;

	/// The counts in each bin.
	std::vector<double> binned_data;

	/**
	 * \brief Calculates the values of the bins (for a particular dimension).
	 *
	 * \param[in] dmin The minimum value (masked) for the dimension.
	 * \param[in] dmax The maximum value (masked) for the dimension.
	 * \param[in] dwidth The width of each bin (in masked coordinates);
	 * \param[in] bstyle The binning style, including the number of bins.
	 * \return The unmasked values for each bin in this dimension.
	 */
	static std::vector<double> bin_values(double dmin, double dmax,
		double dwidth, std::shared_ptr<const BinStyle> bstyle);

public:
	Histogram() = delete;
	
	/**
	 * \brief Constructor specifying the dimensionality of the data.
	 *
	 * \param[in] ndim_ The dimensionality of the data.
	 */
	Histogram(std::size_t ndim_);

	/**
	 * \brief Adds a data element to the histogram.
	 *
	 * \throw std::invalid_argument if the data has the wrong dimensionality.
	 * \throw std::runtime_error if the bins have already been formed.
	 *
	 * \param[in] v The data.
	 */
	void add_data(std::valarray<double> v);

	/**
	 * \brief Bins the data using the specified binning styles for each
	 *    dimension.
	 *
	 * \throw std::invalid_argument if the number of binning styles doesn't
	 *    match the dimensionality of the data.
	 * \throw std::runtime_error if the histogram has already been binned, if
	 *    one dimension doesn't specify a binning style, or if one dimension
	 *    uses a binning style with 0 bins.
	 * \throw std::size_t if a dimenion (its index is thrown) has a null data
	 *    range (all values are the same) and more than one bin is requested.
	 *
	 * \param[in] binstyles The binning styles.
	 */
	void bin_data(
		const std::vector<std::shared_ptr<const BinStyle>> &binstyles);

	/**
	 * \brief Gets an index that iterates over all the bins.
	 *
	 * \throw std::runtime_error if the data has not yet been binned.
	 *
	 * \return The iterator.
	 */
	CounterIndex begin() const;

	/**
	 * \brief Returns the coordinates (the values of the variables in the middle
	 *    of a bin) for the given bin.
	 *
	 * \throw std::runtime_error if the data has not yet been binned.
	 *
	 * \param[in] index The index of the bin.
	 * \return The coordinates of the bin.
	 */
	std::valarray<double> getCoordinates(const CounterIndex &index) const;

	/**
	 * \brief Returns the bin count for the given bin.
	 *
	 * \throw std::runtime_error if the data has not yet been binned.
	 *
	 * \param[in] index The index of the bin.
	 * \return The bin count of the bin.
	 */
	double getBinCount(const CounterIndex &index) const;
};

} // namespace molstat

#endif
