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

// forward declarations
namespace molstat {
class RandomDistribution;
class BinStyle;
}

/**
 * \brief Class that reads the input deck and sets up the molstat::Simulator
 *    object.
 *
 * Preliminary error checking is done in the SimulatorInputParse::readInput()
 * function. This function only processes the lines and stores them for later
 * use. The SimulatorInputParse::create() function actually builds the
 * molstat::Simulator object, performing additional runtime error checking.
 *
 * This separation of responsibility is to provide the user with better error
 * messages, should there be problems.
 */
class SimulatorInputParse
{
private:
	/// Data structure that stores information about models to be created.
	struct ModelInformation
	{
		/// The name of the model to instantiate.
		std::string name;

		/// The list of distributions for this model.
		std::map<std::string,
		          std::shared_ptr<const molstat::RandomDistribution>> dists;

		/// A list of submodels to be created.
		std::list<ModelInformation> submodels;

		/**
		 * \brief Gets a string representation of the model information.
		 *
		 * \return The string.
		 */
		std::string to_string() const;
	};

	/// The top-level simulate model information.
	ModelInformation top_model;

	/**
	 * \brief Map of observable index (axis) to a pair of observable name
	 *    and the binning style.
	 */
	std::map<std::size_t,
	         std::pair<std::string, std::shared_ptr<molstat::BinStyle>>>
		obs_bins;

	/// File name for the histogram output.
	std::string histfilename{ "histogram.dat" };

	/// The number of trials (i.e., data points to simulate).
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
	 * \param[in,out] output Output stream for any error messages.
	 * \param[in,out] lineno The input line number.
	 * \return The model information, save for the name of the model type.
	 */
	static ModelInformation readModel(std::istream &input, std::ostream &output,
		std::size_t &lineno);

	/**
	 * \brief Constructs a model from the
	 *    SimulatorInputParser::ModelInformation.
	 *
	 * \throw std::exception if the model cannot be constructed.
	 *
	 * \param[in,out] output Output stream for any error messages.
	 * \param[in] models Map of available models.
	 * \param[in] info The model information from the input deck.
	 * \return The constructed model.
	 */
	static std::shared_ptr<molstat::SimulateModel> constructModel(
		std::ostream &output,
		const std::map<std::string,
		               molstat::SimulateModelFactoryFunction> &models,
		ModelInformation &info);

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
	 * \param[in,out] output Output stream for any error messages.
	 */
	void readInput(std::istream &input, std::ostream &output);

	/**
	 * \brief Processes the input to actually instantiate models and the
	 *    simulator.
	 *
	 * This function may change the state of the input parser if certain
	 * aspects of the input data are invalid, incomplete, etc.
	 *
	 * \throw exception if any exceptions from the molstat::Simulator or
	 *    molstat::SimulateModel functions is thrown.
	 *
	 * \param[in,out] output Output stream for any error messages.
	 * \return The simulator.
	 */
	std::unique_ptr<molstat::Simulator> createSimulator(std::ostream &output);

	/**
	 * \brief Gets the number of trials.
	 *
	 * \return The number of trials.
	 */
	std::size_t numTrials() const noexcept;

	/**
	 * \brief Prints the state of the input parser.
	 *
	 * For debugging purposes, only.
	 *
	 * \param[in,out] output The output stream.
	 */
	void printState(std::ostream &output) const;

	/**
	 * \brief Returns the name of the output file.
	 *
	 * \return The name of the output file.
	 */
	std::string outputFileName() const;

	/**
	 * \brief Get the binning styles.
	 *
	 * \return A vector containing the binning styles.
	 */
	std::vector<std::shared_ptr<molstat::BinStyle>> getBinStyles() const;
};

#endif
