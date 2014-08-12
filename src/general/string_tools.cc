/* This file is a part of MolStat, which is distributed under the Creative
   Commons Attribution-NonCommercial 4.0 International Public License.
   MolStat (c) 2014, Northwestern University. */
/**
 * \internal
 * \file string_tools.cc
 * \brief Functions for processing strings, useful in I/O.
 *
 * \author Matthew G.\ Reuter
 * \date April 2014
 * \endinternal
 */

#include "string_tools.h"
#include <cstdlib>
#include <stdexcept>
#include <algorithm>

using namespace std;

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

string getline(FILE *f) {
	string line;
	char *cline;
	size_t chars;

	cline = nullptr;
	chars = 0;
	chars = getline(&cline, &chars, f);
	if(chars == -1)
		throw runtime_error("EOF encountered");

	line = cline;
	free(cline);
	return line;
}

void tokenize(const std::string &str, std::vector<std::string> &vec) {
	string::const_iterator next = str.begin();
	string token;

	vec.clear();
	while(next_token(next, str.end(), token))
		vec.push_back(token);
}

void make_lower(std::string &str) {
	for(char &c : str)
		c = tolower(c);
}

double string_to_double(const std::string &str) {
	double ret;

	try {
		ret = stod(str);
	}
	catch(const invalid_argument &e) {
		throw invalid_argument("Unable to convert \"" + str + "\" to a double.");
	}

	return ret;
}
