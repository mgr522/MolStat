/* This file is a part of MolStat, which is distributed under the Creative
   Commons Attribution-NonCommercial 4.0 International Public License.
   MolStat (c) 2014, Northwestern University. */
/**
 * \file string_tools.h
 * \brief Functions for processing strings, useful in I/O.
 *
 * \author Matthew G.\ Reuter
 * \date April 2014
 */

#ifndef __string_tools_h__
#define __string_tools_h__

#include <string>
#include <vector>
#include <cstdio>

/**
 * \internal
 * \brief Gets a line from the desired stream and puts it in std::string form.
 *
 * \exception std::runtime_error if EOF is encountered.
 *
 * \param[in,out] f The stream.
 * \return The string.
 * \endinternal
 */
std::string getline(FILE *f);

/**
 * \internal
 * \brief Tokenizes a string.
 *
 * \param[in] str The string to be tokenize.
 * \param[out] vec The strings, in std::vector form.
 * \endinternal
 */
void tokenize(const std::string &str, std::vector<std::string> &vec);

/**
 * \internal
 * \brief Makes a string lower case.
 *
 * \param[in,out] str The string.
 * \endinternal
 */
void make_lower(std::string &str);

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
 * \brief Turns a string into a double.
 *
 * \exception std::invalid_argument if the string cannot be cast to a double.
 *
 * \param[in] str The string to be cast.
 * \return The double.
 * \endinternal
 */
double string_to_double(const std::string &str);

#endif
