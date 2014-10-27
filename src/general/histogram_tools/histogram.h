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
 *
 * An iterator-style class is provided for iterating through the bins.
 */
class Histogram {
#if 0
public:
	/**
	 * \internal
	 * \brief Iterator class for accessing the histogram bins.
	 *
	 * This class needs the inverse mask function for calculating the middle of
	 * the bin.
	 * \endinternal
	 */
	class const_iterator {
	protected:
		/**
		 * \internal
		 * \brief The bin index.
		 * \endinternal
		 */
		std::valarray<size_t> bin;

		/**
		 * \internal
		 * \brief The value of the variable(s) in the middle of this bin.
		 * \endinternal
		 */
		std::valarray<double> val;

		/**
		 * \internal
		 * \brief The bin count of this bin.
		 * \endinternal
		 */
		double bincount;

		/**
		 * \internal
		 * \brief The binning styles.
		 * \endinternal
		 */
		std::vector<std::shared_ptr<const BinStyle>> bstyles;

		/**
		 * \internal
		 * \brief Advances the values of bin to the next bin.
		 * \endinternal
		 */
		void next_bin();

		/**
		 * \internal
		 * \brief Sets the output values for this bin.
		 * \endinternal
		 */
		void set_output();

	public:
		const_iterator() = delete;

		/**
		 * \internal
		 * \brief Constructor specifying the binning style for each axis.
		 *
		 * \param[in] bstyle_ The binning style.
		 * \endinternal
		 */
		const_iterator(const std::shared_ptr<const BinStyle> bstyle_);

		/**
		 * \internal
		 * \brief Get the value of the variable in the middle of this bin.
		 *
		 * \return The value of the variable in the middle of this bin.
		 * \endinternal
		 */
		const std::valarray<double> &get_variable() const;

		/**
		 * \internal
		 * \brief Get the bin count of this bin.
		 *
		 * \return The bin count of this bin.
		 * \endinternal
		 */
		double get_bin_count() const;
	};
#endif

private:
	/**
	 * \brief State of the histogram.
	 *
	 * False if binning has not occurred; true if it has.
	 */
	bool haveBinned;

	/**
	 * \brief The dimensionality of the data.
	 */
	const std::size_t ndim;

	/**
	 * \brief The accumulated data.
	 */
	std::forward_list<std::valarray<double>> data;

	/**
	 * \brief The minimum and maximum values in each dimension.
	 */
	std::vector<std::array<double, 2>> extremes;

	/**
	 * \brief The middle of each bin.
	 *
	 * The first index is the dimension, the second index is the bin number.
	 */
	std::vector<std::vector<double>> bin_value;

	/**
	 * \brief The counts in each bin.
	 */
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
	 * \throw std::runtime_error if one dimension does not have a range of
	 *    values and more than one bin is requested.
	 *
	 * \param[in] binstyles The binning styles.
	 */
	void bin_data(
		const std::vector<std::shared_ptr<const BinStyle>> &binstyles);
};

} // namespace molstat

#endif
