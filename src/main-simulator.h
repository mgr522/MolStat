/* This file is a part of MolStat, which is distributed under the Creative
   Commons Attribution-NonCommercial 4.0 International Public License.
   MolStat (c) 2014, Northwestern University. */
/**
 * \file main-simulator.h
 * \brief Function and class declarations for the main simulator program.
 *
 * \author Matthew G.\ Reuter
 * \date October 2014
 */

#ifndef __main_simulator_h__
#define __main_simulator_h__

#include <iostream>
#include <memory>
#include <string>
#include <queue>
#include <map>

#include <general/simulator_tools/simulator.h>

/**
 * \internal
 * \brief Enum for the type of histogram (1D or 2D).
 * \endinternal
 */
enum class HistogramType
{
	OneD, // 1D Histogram
	TwoD // 2D Histogram
};

/**
 * \brief Processes the input deck.
 *
 * \throw std::runtime_error if an unrecoverable error occurs.
 *
 * \param[in,out] input The input stream.
 * \param[out] ofname Name for the histogram output file.
 * \return Pointer to the molstat::Simulator object.
 */
std::unique_ptr<molstat::Simulator> processInput(std::istream &input,
	std::string &ofname);

/**
 * \brief Gets a molstat::SimulateModel from the tokens.
 *
 * The command for making a model is multi-line; this function has input loop
 * similar to that of ::processInput.
 *
 * The list of tokens is destroyed in this function.
 *
 * \throw std::runtime_error if
 *    - no model name is specified,
 *    - the model name is not found,
 *    - EOF is encountered while reading the model's command block.
 *
 * \param[in] tokens The list of tokens for the first line of the model command.
 * \param[in] input The input stream.
 * \param[in] models The list of available models (stored as a map from string
 *    to a function that facilitates model creation).
 * \return The constructed model.
 */
std::shared_ptr<molstat::SimulateModel> processModel(
	std::queue<std::string> &&tokens,
	std::istream &input,
	const std::map<std::string,
	               molstat::SimulateModelFactoryFunction> &models);

/**
 * \brief Gets an observable from a list of tokens.
 *
 * The list of tokens is destroyed in this function.
 *
 * \throw std::runtime_error if the token cannot be matched to an observable or
 *    a name was not specified (too few tokens).
 *
 * \param[in] tokens The list of tokens.
 * \param[in] observables The list of observables (stored as a map from
 *    string to ObservableIndex).
 * \return The observable.
 */
molstat::ObservableIndex processObservable(
	std::queue<std::string> &&tokens,
	const std::map<std::string,
	               molstat::ObservableIndex> &observables);

#endif