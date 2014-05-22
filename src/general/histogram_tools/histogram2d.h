/* This file is a part of MolStat, which is distributed under the Creative
   Commons Attribution-NonCommercial 4.0 International Public License.
   MolStat (c) 2014, Northwestern University. */
/**
 * \file histogram2d.h
 * \brief Provides a class for making 2D histograms from data.
 *
 * \author Matthew G.\ Reuter
 * \date May 2014
 */

#ifndef __histogram2d_h__
#define __histogram2d_h__

#include <memory>
#include <gsl/gsl_histogram2d.h>
#include "histogram_interface.h"
#include "bin_style.h"

/**
 * \brief Implements 2D histograms.
 *
 * An iterator-style class is also provided for iterating through the bins.
 */
class Histogram2D : public Histogram<2> {
protected:
	/**
	 * \brief The handle to the GSL histogram functions.
	 */
	std::shared_ptr<gsl_histogram2d> hist;

public:
	class const_iterator;
	Histogram2D() = delete;

	/**
	 * \brief Constructor requiring the number of bins in the histogram, the
	 *    ranges of the variables, and the binning style.
	 *
	 * \param[in] nbin The number of bins to use, in each dimension.
	 * \param[in] mins The minimum values in the ranges of each dimension.
	 * \param[in] maxs The maximum values in the ranges of each dimension.
	 * \param[in] bstyle The binning style.
	 */
	Histogram2D(const std::array<std::size_t, 2> &nbin,
		const std::array<double, 2> &mins, const std::array<double, 2> &maxs,
		const std::shared_ptr<const BinStyle> bstyle);

	/**
	 * \brief Adds a data element to the histogram.
	 *
	 * \param[in] v The 2D data variable.
	 */
	virtual void add_data(const std::array<double, 2> &v);

	/**
	 * \brief Get an iterator to the front of this histogram.
	 *
	 * \return An iterator to the beginning of this histogram.
	 */
	const_iterator begin() const;

	/**
	 * \brief Get an iterator to the end of this histogram.
	 *
	 * \return An iterator to the end of this histogram.
	 */
	const_iterator end() const;

	/**
	 * \brief Iterator class for accessing the histogram bins.
	 */
	class const_iterator : public Histogram<2>::const_iterator {
	protected:
		/**
		 * \brief The GSL histogram handle.
		 */
		const std::shared_ptr<const gsl_histogram2d> hist;

		/**
		 * \brief Advances the values of bin1 and bin2 to the next bin.
		 */
		virtual void next_bin();

		/**
		 * \brief Sets the values of val1_, val2_, and bincount_ for this bin.
		 */
		virtual void set_output();

	public:
		const_iterator() = delete;

		/**
		 * \brief Constucts the iterator; based on the GSL histogram handle.
		 *
		 * \param[in] h The GSL histogram handle.
		 * \param[in] bstyle The binning style.
	 	 */
		const_iterator(const std::shared_ptr<const gsl_histogram2d> h,
			const std::shared_ptr<const BinStyle> bstyle);

		/**
		 * \brief Prefix forward iteration operator.
		 *
		 * \return State of the iterator before iteration.
		 */
		const_iterator operator++();

		/**
		 * \brief Postfix forward iteration operator.
		 *
		 * \return State of the iterator after iteration.
		 */
		const_iterator operator++(int);

		/**
		 * \brief Determine if these iterators refer to the same state.
		 *
		 * \param[in] rhs The iterator to compare to.
		 * \return True if they have the same state; false otherwise.
		 */
		bool operator== (const const_iterator &rhs) const;

		/**
		 * \brief Determine if these iterators refer to different states.
		 *
		 * \param[in] rhs The iterator to compare to.
		 * \return True if they have different states; false otherwise.
		 */
		bool operator!= (const const_iterator &rhs) const;

		/**
		 * \brief Histogram2D::end() const needs to set the iterator to the end.
	 	 */
		friend const_iterator Histogram2D::end() const;
	};
};

#endif
