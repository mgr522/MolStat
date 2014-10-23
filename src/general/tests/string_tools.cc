/* This file is a part of MolStat, which is distributed under the Creative
   Commons Attribution-NonCommercial 4.0 International Public License.
   MolStat (c) 2014, Northwestern University. */
/**
 * \internal
 * \file tests/string_tools.cc
 * \brief Test suite for the string functions.
 *
 * \test Tests of the various string handling functions.
 *
 * \author Matthew G.\ Reuter
 * \date October 2014
 * \endinternal
 */

#include <list>
#include <queue>
#include <cassert>
#include <string>
#include <general/string_tools.h>

using namespace std;

/**
 * \internal
 * \brief Auxiliary function that makes sure the tokens are the same.
 *
 * The containers are destroyed on exit. q2 is a list because a queue
 * cannot be constructed from an initializer_list.
 *
 * \param[in] q1 A queue of tokens (presumably from molstat::tokenize).
 * \param[in] q2 A list of tokens to compare against.
 * \endinternal
 */
static void assert_tokens_same(queue<string> &&q1,
	list<string> &&q2)
{

	assert(q1.size() == q2.size());
	while(q1.size() > 0)
	{
		assert(q1.front() == q2.front());
		q1.pop();
		q2.pop_front();
	}
}

/**
 * \internal
 * \brief Main function for testing the various string functions.
 *
 * \param[in] argc The number of command-line arguments.
 * \param[in] argv The command-line arguments.
 * \return Exit status: 0 if the code passes the test, non-zero otherwise.
 * \endinternal
 */
int main(int argc, char **argv)
{
	// test the to_lower function
	assert(string("") == molstat::to_lower(""));
	assert(string("hello") == molstat::to_lower("hello"));
	assert(string("hello") == molstat::to_lower("HeLLo"));
	assert(string("Oops") != molstat::to_lower("Oops"));
	assert(string("12345asdfg") == molstat::to_lower("12345ASDfg"));

	// test the tokenize function
	assert_tokens_same(
		molstat::tokenize("Hello, world!"),
		{ "Hello,", "world!" });
	
	assert_tokens_same(
		molstat::tokenize("A string with a \"multi-word phrase\" in it."),
		{ "A", "string", "with", "a", "multi-word phrase", "in", "it." });

	assert_tokens_same(
		molstat::tokenize("    \tLeading whitespace"),
		{ "Leading", "whitespace" });

	assert_tokens_same(
		molstat::tokenize(" Other forms\tof\nwhitespace   "),
		{ "Other", "forms", "of", "whitespace" });

	// check the cast function
	// to size_t
	assert(4 == molstat::cast_string<size_t>("4"));
	assert(0 == molstat::cast_string<size_t>("0"));

	// some bad size_t casts
	try
	{
		molstat::cast_string<size_t>("-1");

		assert(false);
	}
	catch(const bad_cast &e)
	{
		// should be here
	}

	try
	{
		molstat::cast_string<size_t>("a");
		assert(false);
	}
	catch(const bad_cast &e)
	{
		// should be here
	}

	try
	{
		molstat::cast_string<size_t>("5-1");
		assert(false);
	}
	catch(const bad_cast &e)
	{
		// should be here
	}

	// to double
	assert(abs(4.5 - molstat::cast_string<double>("4.5")) < 1.e-6);
	assert(abs(1.1e2 - molstat::cast_string<double>("1.1e2")) < 1.e-6);
	assert(abs(1.1e2 - molstat::cast_string<double>("1.1E2")) < 1.e-6);
	assert(abs(-4. - molstat::cast_string<double>("-4.000000")) < 1.e-6);
	assert(abs(-4. - molstat::cast_string<double>("-4")) < 1.e-6);

	// some bad double casts
	try
	{
		molstat::cast_string<double>("a");	
		assert(false);
	}
	catch(const bad_cast &e)
	{
		// should be here
	}

	try
	{
		molstat::cast_string<double>("_-");
		assert(false);
	}
	catch(const bad_cast &e)
	{
		// should be here
	}

	try
	{
		// this one fails because the entire string/token is not used
		// std::stod would simply return 4.5
		molstat::cast_string<double>("4.5-1.4");
		assert(false);
	}
	catch(const bad_cast &e)
	{
		// should be here
	}

	// unimplemented cast
	try
	{
		 molstat::cast_string<queue<list<string>>>("abcdefg");
		 assert(false);
	}
	catch(const bad_cast &e)
	{
		// should be here
	}
	
	return 0;
}
