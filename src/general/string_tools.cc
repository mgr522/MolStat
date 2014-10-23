/* This file is a part of MolStat, which is distributed under the Creative
   Commons Attribution-NonCommercial 4.0 International Public License.
   MolStat (c) 2014, Northwestern University. */
/**
 * \internal
 * \file string_tools.cc
 * \brief Functions for processing strings, useful in I/O.
 *
 * \author Matthew G.\ Reuter
 * \date October 2014
 * \endinternal
 */

#include "string_tools.h"
#include <cstdlib>
#include <stdexcept>
#include <algorithm>

using namespace std;

namespace molstat {

/**
 * \internal
 * \brief Gets the next token in the string.
 *
 * \param[in,out] next The next character to examine (via iterator).
 * \param[in] end The end of the string.
 * \param[out] tok The token, if a valid token exists.
 * \return True if a token is found, false otherwise.
 * \endinternal
 */
static bool next_token(string::const_iterator &next,
	string::const_iterator end, string &tok) {

	// skip past any leading whitespace
	while(next != end && isspace(*next))
		++next;
	if(next == end)
		return false;

	if(*next == '"') {
		// this token is delimited by quotes
		string::const_iterator quote = next;
		do {
			++quote; // move past the current quote mark
			quote = find(quote, end, '"');
			if(quote == end) // unterminated token
				return false;
		} while(*(quote-1) == '\\'); // allow for escaped quotes inside

		tok.assign(next+1, quote); // this omits the quotes
		next = quote;
		++next; // move past the terminal quote
	}
	else if(*next == '<') {
		// this token is delimited by angle brackets
		string::const_iterator bracket = find(next, end, '>');
		if(bracket == end) // unterminated token
			return false;
		tok.assign(next, bracket+1);
		next = bracket+1; // move past the closing bracket
	}
	else {
		string::const_iterator first = next;
		while(next != end && !isspace(*next))
			++next;
		tok.assign(first, next);
	}

	return true;
}

TokenContainer tokenize(const std::string &str) {
	TokenContainer ret;
	string::const_iterator next = str.begin();
	string token;

	while(next_token(next, str.end(), token))
		ret.emplace(token);

	return ret;
}

std::string to_lower(const std::string &str) {
	std::string ret{ str };

	for(char &c : ret)
		c = tolower(c);

	return ret;
}

template<>
double cast_string(const std::string &str) {
	double ret;
	size_t next;

	try {
		// do the conversion and get the index of the first character not used.
		ret = stod(str, &next);

		// make sure the entire string was used.
		if(next != str.size())
			throw bad_cast();
	}
	catch(const invalid_argument &e) {
		throw bad_cast();
	}

	return ret;
}

} // namespace molstat