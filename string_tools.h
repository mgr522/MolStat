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

/**
 * \brief Tokenizes a string.
 *
 * \param[in] str The string to be tokenize.
 * \param[out] vec The strings, in std::vector form.
 */
void tokenize(const std::string &str, std::vector<std::string> &vec);

/**
 * \brief Makes a string lower case.
 *
 * \param[in,out] str The string.
 */
void make_lower(std::string &str);

#endif
