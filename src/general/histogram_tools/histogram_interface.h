/* This file is a part of MolStat, which is distributed under the Creative
   Commons Attribution-NonCommercial 4.0 International Public License.
   MolStat (c) 2014, Northwestern University. */
/**
 * \file histogram_interface.h
 * \brief Provides an interface for the GSL histogram (1D and 2D) functions.
 *
 * \author Matthew G.\ Reuter
 * \date May 2014
 */

#ifndef __histogram_interface_h__
#define __histogram_interface_h__

#include <memory>
#include <array>
#include "bin_style.h"

/**
 * \internal
 * \brief Interface for making histograms using GSL.
 *
 * This class provides a consistent interface to the GSL histogram functions
 * (both 1D and 2D). At construction, the user specifies the number of bins in
 * the histogram, as well as the minimum and maximum values for all variables.
 * The data elements are then added (one at a time) to make the histogram.
 *
 * The user also needs to specify the binning style for how the data is to be
 * binned. More information on binning styles can be found in the notes for
 * BinningStyle.
 *
 * An iterator-style class is provided for iterating through the bins.
 *
 * \tparam N The dimensionality of the histogram.
 * \endinternal
 */
template <std::size_t N>
class Histogram {
protected:
	/**
	 * \internal
	 * \brief Iterator class for accessing the histogram bins.
	 *
	 * This class needs the inverse conductance function for calculating the
	 * middle of the bin.
	 * \endinternal
	 */
	class const_iterator {
	protected:
		/**
		 * \internal
		 * \brief The bin index.
		 * \endinternal
		 */
		std::array<size_t, N> bin;

		/**
		 * \internal
		 * \brief The value of the variable(s) in the middle of this bin.
		 * \endinternal
		 */
		std::array<double, N> val;

		/**
		 * \internal
		 * \brief The bin count of this bin.
		 * \endinternal
		 */
		double bincount;

		/**
		 * \internal
		 * \brief The binning style.
		 * \endinternal
		 */
		std::shared_ptr<const BinStyle> bstyle;

		/**
		 * \internal
		 * \brief Advances the values of bin to the next bin.
		 * \endinternal
		 */
		virtual void next_bin() = 0;

		/**
		 * \internal
		 * \brief Sets the output values for this bin.
		 * \endinternal
		 */
		virtual void set_output() = 0;

	public:
		const_iterator() = delete;

		/**
		 * \internal
		 * \brief Constructor specifying binning style.
		 *
		 * \param[in] bstyle_ The binning style.
		 * \endinternal
		 */
		const_iterator(const std::shared_ptr<const BinStyle> bstyle_);

		/**
		 * \internal
		 * \brief Default destructor.
		 * \endinternal
		 */
		virtual ~const_iterator() = default;

		/**
		 * \internal
		 * \brief Get the value of the variable in the middle of this bin.
		 *
		 * \return The value of the variable in the middle of this bin.
		 * \endinternal
		 */
		const std::array<double, N> &get_variable() const;

		/**
		 * \internal
		 * \brief Get the bin count of this bin.
		 *
		 * \return The bin count of this bin.
		 * \endinternal
		 */
		double get_bin_count() const;
	};

	/**
	 * \internal
	 * \brief The binning style.
	 * \endinternal
	 */
	std::shared_ptr<const BinStyle> bstyle;

public:
	Histogram() = delete;

	/**
	 * \internal
	 * \brief Constructor that processes the binning style.
	 *
	 * \param[in] bstyle_ The binning style.
	 * \endinternal
	 */
	Histogram(const std::shared_ptr<const BinStyle> bstyle_);

	/**
	 * \internal
	 * \brief Destructor.
	 * \endinternal
	 */
	virtual ~Histogram() = default;

	/**
	 * \internal
	 * \brief Adds a data element to the histogram.
	 *
	 * \param[in] v The variable.
	 * \endinternal
	 */
	virtual void add_data(const std::array<double, N> &v) = 0;
};

// Implementations of non-virtual functions
template<std::size_t N>
Histogram<N>::Histogram(const std::shared_ptr<const BinStyle> bstyle_)
	: bstyle(bstyle_) {
}

template<std::size_t N>
Histogram<N>::const_iterator::const_iterator(
	const std::shared_ptr<const BinStyle> bstyle_)
	: bincount(0), bstyle(bstyle_) {

	for(std::size_t i = 0; i < N; ++i) {
		bin[i] = 0;
		val[i] = 0.;
	}
}

template<std::size_t N>
const std::array<double, N> &Histogram<N>::const_iterator::get_variable() const
{
	return val;
}

template<std::size_t N>
double Histogram<N>::const_iterator::get_bin_count() const {
	return bincount;
}

#endif
