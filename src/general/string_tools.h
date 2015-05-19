/* This file is a part of MolStat, which is distributed under the Creative
   Commons Attribution-NonCommercial 4.0 International Public License.

   (c) 2014 Northwestern University. */

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
#include <typeinfo>

/**
 * \namespace molstat
 * \brief Namespace for components of MolStat.
 */
namespace molstat {

/// Alias for a container of tokens.
using TokenContainer = std::queue<std::string>;

/**
 * \brief Tokenizes a string.
 *
 * Tokens are delimited by whitespace. Quotes (\"\") can be used to include
 * whitespace in a token.
 *
 * \param[in] str The string to be tokenized.
 * \return The tokens.
 */
TokenContainer tokenize(const std::string &str);

/**
 * \brief Returns a copy of the string, in lower case.
 *
 * \param[in] str The string.
 * \return The string, in lower case.
 */
std::string to_lower(const std::string &str);

/**
 * \brief Finds and replaces a substring within another string.
 *
 * \param[in] str String modified by the search and replace.
 * \param[in] find The pattern to replace.
 * \param[in] replace The string to substitute.
 * \return The substitued string.
 */
std::string find_replace(const std::string &str, const std::string &find,
	const std::string &replace);

/**
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
 */
template<typename T>
T cast_string(const std::string &str) {
	throw std::bad_cast();
}

/**
 * \brief Specialization for casting a string to a double.
 *
 * The entire string must be used in the conversion.
 *
 * \throw std::bad_cast if the string cannot be cast to the double or if only
 *    the beginning of the string is used.
 *
 * \param[in] str The string to be cast.
 * \return The string, as a double.
 */
template<>
double cast_string(const std::string &str);

/**
 * \brief Specialization for casting a string to a size_t.
 *
 * The entire string must be used in the conversion.
 *
 * \throw std::bad_cast if the string cannot be cast to the double or if only
 *    the beginning of the string is used.
 *
 * \param[in] str The string to be cast.
 * \return The string, as a double.
 */
template<>
std::size_t cast_string(const std::string &str);

} // namespace molstat

#endif
