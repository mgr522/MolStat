/**
 * \file histogram1d.h
 * \brief Provides an interface for the GSL histogram functions.
 *
 * \author Matthew G.\ Reuter
 * \date May 2014
 */

#ifndef __histogram1d_h__
#define __histogram1d_h__

#include <memory>
#include <gsl/gsl_histogram.h>
#include "histogram_interface.h"

/**
 * \brief Implements 1D histograms.
 *
 * An iterator-style class is also provided for iterating through the bins.
 */
class Histogram1D : public Histogram<1> {
protected:
	/**
	 * \brief The handle to the GSL histogram functions.
	 */
	std::shared_ptr<gsl_histogram> hist;

public:
	class const_iterator;

	/**
	 * \brief Constructor requiring the number of bins in the histogram and the
	 *    ranges of the variables.
	 *
	 * \param[jn] nbin The number of bins to use.
	 * \param[in] minval The minimum value in the histogram range.
	 * \param[in] maxval The maximum value in the histogram range.
	 */
	Histogram1D(const std::size_t nbin, const double minval,
		const double maxval);

	/**
	 * \brief Adds a data element to the histogram.
	 *
	 * \param[in] v The variable.
	 */
	void add_data(const double v);

	/**
	 * \brief Adds a data element to the histogram.
	 *
	 * \param[in] v The variable.
	 */
	virtual void add_data(const std::array<double, 1> &v);

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
	class const_iterator : public Histogram<1>::const_iterator {
	protected:
		/**
		 * \brief The GSL histogram handle.
		 */
		const std::shared_ptr<const gsl_histogram> hist;

		/**
		 * \brief Advances the values of bin to the next bin.
		 */
		virtual void next_bin();

		/**
		 * \brief Sets the values of val_ and bincount_ for this bin.
		 */
		virtual void set_output();

	public:
		const_iterator() = delete;

		/**
		 * \brief Constucts the iterator; based on the GSL histogram handle.
		 *
		 * \param[in] h The GSL histogram handle.
	 	 */
		const_iterator(const std::shared_ptr<const gsl_histogram> h);

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
		 * \brief Histogram1D::end() const needs to set the iterator to the end.
	 	 */
		friend const_iterator Histogram1D::end() const;
	};
};

#endif
