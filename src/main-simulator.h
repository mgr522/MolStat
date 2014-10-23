/* This file is a part of MolStat, which is distributed under the Creative
   Commons Attribution-NonCommercial 4.0 International Public License.
   MolStat (c) 2014, Northwestern University. */
/**
 * \file main-simulator.h
 * \brief Function and class declarations for the main simulator program and
 *    the molstat::SimulatorInputParse class for parsing the input deck.
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
#include <list>
#include <map>

#include <general/simulator_tools/simulator.h>

// forward declaration
namespace molstat {
class RandomDistribution;
class BinStyle;
}

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
 * \brief Class that reads the input deck and sets up the molstat::Simulator
 *    object.
 *
 * Preliminary error checking is done in the SimulatorInputParse::readInput()
 * function. This function only process the lines and stores them for later
 * use. The SimulatorInputParse::create() function actually builds the
 * molstat::Simulator object, performing additional runtime error checking.
 *
 * This separation of responsibility is to provide the user with better error
 * messages, should there be problems.
 */
class SimulatorInputParse {
private:
	/**
	 * \brief Data structure that stores information about models to be created.
	 */
	struct ModelInformation {
		/**
		 * \brief The name of the model to instantiate.
		 */
		std::string name;

		/**
		 * \brief The list of distributions for this model.
		 */
		std::map<std::string,
		          std::shared_ptr<const molstat::RandomDistribution>> dists;

		/**
		 * \brief A list of submodels to be created.
		 */
		std::list<ModelInformation> submodels;

		/**
		 * \internal
		 * \brief Gets a string representation of the model information.
		 *
		 * \return The string.
		 * \endinternal
		 */
		std::string to_string() const;
	};

	/**
	 * \brief The top-level simulate model information.
	 */
	ModelInformation top_model;

	/**
	 * \brief Map of observable index (axis) to the observable name.
	 */
	std::map<std::size_t, std::string> used_observables;

	/**
	 * \brief Map of observable index (axis) to the binning style tokens.
	 */
	std::map<std::size_t, std::shared_ptr<molstat::BinStyle>> bin_styles;

	/**
	 * \brief File name for the histogram output.
	 */
	std::string histfilename{ "histogram.dat" };

	/**
	 * \brief The number of trials (i.e., data points to simulate).
	 */
	std::size_t trials{ 0 };

	/**
	 * \brief Prints an error message.
	 *
	 * \param[in,out] output The output stream.
	 * \param[in] lineno The line number.
	 * \param[in] message The error message.
	 */
	static inline void printError(std::ostream &output, std::size_t lineno,
		std::string message);

	/**
	 * \brief Reads a model from the input stream
	 *
	 * \throw std::runtime_error if there is a fatal error reading the input
	 *    deck (notably, a lack of "endmodel" before EOF).
	 *
	 * \param[in,out] input The input stream.
	 * \param[in,out] lineno The input line number.
	 * \return The model information, save for the name of the model type.
	 */
	static ModelInformation readModel(std::istream &input, std::size_t &lineno);

	/**
	 * \brief Constructs a model from the
	 *    SimulatorInputParser::ModelInformation.
	 *
	 * \throw std::exception if the model cannot be constructed.
	 *
	 * \param[in] models Map of available models.
	 * \param[in] info The model information from the input deck.
	 * \return The constructed model.
	 */
	static std::shared_ptr<molstat::SimulateModel> constructModel(
		const std::map<std::string,
		               molstat::SimulateModelFactoryFunction> &models,
		const ModelInformation &info);

public:
	/**
	 * \brief Reads the input deck from the stream and performs some runtime
	 *    error checking.
	 *
	 * This function minimally processes the input deck; it does not attempt
	 * to instantiate models or the simulator. Instead, it creates a data
	 * structure that can be processed in the next call.
	 *
	 * \throw std::runtime_error if there is a fatal error reading the input
	 *    deck.
	 *
	 * \param[in,out] input The input stream containing the input deck.
	 */
	void readInput(std::istream &input);

	/**
	 * \brief Processes the input to actually instantiate models and the
	 *    simulator.
	 *
	 * \throw exception if any exceptions from the molstat::Simulator or
	 *    molstat::SimulateModel functions is thrown.
	 *
	 * \return The simulator.
	 */
	std::unique_ptr<molstat::Simulator> createSimulator();

	/**
	 * \brief Gets the number of trials.
	 *
	 * \return The number of trials.
	 */
	std::size_t numTrials() const noexcept;

	/**
	 * \internal
	 * \brief Prints the state of the input parser.
	 *
	 * For debugging purposes, only.
	 *
	 * \param[in,out] output The output stream.
	 * \endinteral
	 */
	void printState(std::ostream &output) const;
};

#if 0
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

#endif