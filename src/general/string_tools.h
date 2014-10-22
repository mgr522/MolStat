/* This file is a part of MolStat, which is distributed under the Creative
   Commons Attribution-NonCommercial 4.0 International Public License.
   MolStat (c) 2014, Northwestern University. */
/**
 * \file string_tools.h
 * \brief Functions for processing strings, useful in I/O.
 *
 * \author Matthew G.\ Reuter
 * \date October 2014
 */

#ifndef __string_tools_h__
#define __string_tools_h__

#include <string>
#include <queue>

/**
 * \namespace molstat
 * \brief Namespace for components of MolStat.
 */
namespace molstat {

/**
 * \internal
 * \brief Tokenizes a string.
 *
 * Tokens are delimited by whitespace. Quotes (\"\") can be used to include
 * whitespace in a token.
 *
 * \param[in] str The string to be tokenized.
 * \return The tokens, ordered in a std::queue.
 * \endinternal
 */
std::queue<std::string> tokenize(const std::string &str);

/**
 * \internal
 * \brief Returns a copy of the string, in lower case.
 *
 * \param[in] str The string.
 * \return The string, in lower case.
 */
std::string to_lower(const std::string &str);

/**
 * \internal
 * \brief Template for casting a string to the desired type.
 *
 * The entire string must be used in the conversion.
 *
 * \throw std::bad_cast if the string cannot be cast to a specified type or if
 *    only part of the string was needed.
 *
 * \tparam T The type to which we aim to cast.
 * \param[in] str The string to be cast.
 * \return The string, as cast to the specified type.
 * \endinternal
 */
template<typename T>
T cast_string(const std::string &str) {
	throw std::bad_cast();
}

/**
 * \internal
 * \brief Specialization for casting a string to a double.
 *
 * The entire string must be used in the conversion.
 *
 * \throw std::bad_cast if the string cannot be cast to the double or if only
 *    the beginning of the string is used.
 *
 * \param[in] str The string to be cast.
 * \return The string, as a double.
 * \endinternal
 */
template<>
double cast_string(const std::string &str);

} // namespace molstat

#endif
