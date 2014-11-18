/* This file is a part of MolStat, which is distributed under the Creative
   Commons Attribution-NonCommercial 4.0 International Public License.
   MolStat (c) 2014, Northwestern University. */
/**
 * \file main-simulator-inputparse.cc
 * \brief Definition of the SimulatorInputParse class for parsing the simulator
 *    input deck.
 *
 * \todo Document how this code works.
 *
 * \author Matthew G.\ Reuter
 * \date October 2014
 */

#include "main-simulator.h"
#include <iomanip>
#include <general/simulator_tools/simulator_exceptions.h>
#include <general/random_distributions/rng.h>
#include <general/histogram_tools/bin_style.h>
#include <general/simulator_tools/identity_tools.h>
#include <electron_transport/simulator_models/transport_simulate_module.h>

using namespace std;

inline void SimulatorInputParse::printError(std::ostream &output,
	std::size_t lineno, std::string message)
{
	output << "Error on line " << setw(2) << lineno << ": " << message << endl;
}

std::unique_ptr<molstat::Simulator> SimulatorInputParse::createSimulator(
	std::ostream &output)
{
	// the map of model instantiators
	map<string, molstat::SimulateModelFactoryFunction> models;

	// the map of observable indexes
	map<string, molstat::ObservableIndex> observables;

	// load model and observable names
	// the general syntax for loading models is
	// models.emplace( molstat::to_lower(name),
	//                 molstat::GetSimulateModelFactory<model_type>() );
	//
	// and likewise for observables,
	// observables.emplace( molstat::to_lower(name),
	//                      molstat::GetObservableIndex<observable_type>() );
	models.emplace( molstat::to_lower("IdentityModel"),
		molstat::GetSimulateModelFactory<molstat::IdentityModel>() );
	observables.emplace( molstat::to_lower("Identity"),
		molstat::GetObservableIndex<molstat::IdentityObservable>() );

	// load transport models
	molstat::transport::load_models(models);
	molstat::transport::load_observables(observables);

	// make the model
	// if there are exceptions, let them pass up to the caller
	shared_ptr<molstat::SimulateModel> model
		{ constructModel(output, models, top_model) };
	
	// make the simulator
	unique_ptr<molstat::Simulator> sim{ new molstat::Simulator(model) };

	// set the observables
	auto obs_iter = obs_bins.begin();
	while(obs_iter != obs_bins.end())
	{
		// flag for a successful set of the observable
		bool good_set{ false };

		try
		{
			// get the index
			auto obsindex = observables.at(obs_iter->second.first);
			try
			{
				// set the index
				sim->setObservable(obs_iter->first, obsindex);
				good_set = true;
			}
			catch(const exception &e) // problem setting the observable
			{
				output << "Error setting observable " << obs_iter->first <<
					".\n   " << e.what() << endl;
			}
		}
		catch(const out_of_range &e) // index not found
		{
			output << "Unknown observable: \"" << obs_iter->second.first << "\"."
				<< endl;
		}

		if(good_set)
			++obs_iter;
		else // couldn't set the observable
		{
			// remove the observable from the list
			auto obs_here = obs_iter;
			++obs_iter;
			obs_bins.erase(obs_here);
		}
	}

	return sim;
}

void SimulatorInputParse::readInput(std::istream &input, std::ostream &output)
{
	std::size_t lineno{ 1 }; // line number

	// process the input deck
	// getline evaluates to "false" when we hit EOF
	for(string line; getline(input, line);)
	{
		// tokenize the string
		molstat::TokenContainer tokens = molstat::tokenize(line);
		if(tokens.size() == 0) // empty line
		{
			++lineno;
			continue;
		}

		// the first token is the command name, pop it off and then process the
		// rest of the tokens
		string command { molstat::to_lower(tokens.front()) };
		tokens.pop();

		// go through the supported commands
		if(command == "model")
		{
			// the total lineno needs to be passed to the model reader
			const size_t mylineno{ lineno };

			// enter the model reader to process input lines until the
			// appropriate "endmodel" command is found
			ModelInformation model = readModel(input, output, ++lineno);

			// make sure there is a name (model type) specified
			if(tokens.size() == 0)
			{
				printError(output, mylineno, "No model type specified.");
			}
			else
			{
				// store the name of the model type
				model.name = molstat::to_lower(tokens.front());

				// set this to the parser's top model
				top_model = move(model);
			}
		}
		else if(command == "observable" || command == "observable_x" ||
			command == "observable_y")
		{
			// make sure we have a name for the observable, the number of bins,
			// and (at least) the name of the binning style
			if(tokens.size() < 3)
			{
				printError(output, lineno, "No observable, number of bins, " \
					"and/or binning style specified.");
			}
			else
			{
				// store the name of the observable for later
				string obsname = molstat::to_lower(tokens.front());
				tokens.pop();

				// construct the binning style
				try
				{
					shared_ptr<molstat::BinStyle> binstyle
						{ molstat::BinStyleFactory(move(tokens)) };

					// store the observable name and binning style
					if(command == "observable_y")
						obs_bins.emplace(1, make_pair(obsname, binstyle));
					else
						obs_bins.emplace(0, make_pair(obsname, binstyle));
				}
				catch(const invalid_argument &e)
				{
					// indent the error message
					printError(output, lineno,
						molstat::find_replace(e.what(), "\n", "\n   "));
				}
			}
		}
		else if(command == "output")
		{
			if(tokens.size() == 0)
			{
				printError(output, lineno, "No output file name specified.");
			}
			else
			{
				histfilename = tokens.front();
			}
		}
		else if(command == "trials")
		{
			if(tokens.size() == 0)
			{
				printError(output, lineno, "Number of trials not specified.");
			}
			else
			{
				try
				{
					trials = molstat::cast_string<size_t>(tokens.front());
					if(trials == 0)
						printError(output, lineno,
							"More than 0 trials should be specified.");
				}
				catch(const bad_cast &e)
				{
					printError(output, lineno, "Unable to convert \"" + tokens.front() +
						"\" to a non-negative number.");
				}
			}
		}
		else
		{
			printError(output, lineno, "Unknown command: \"" + command + "\".");
		}

		// move to the next line
		++lineno;
	}
}

SimulatorInputParse::ModelInformation SimulatorInputParse::readModel(
	std::istream &input, std::ostream &output, std::size_t &lineno)
{
	ModelInformation ret;

	while(input)
	{
		string line;
		getline(input, line);

		// tokenize the string
		molstat::TokenContainer tokens = molstat::tokenize(line);
		if(tokens.size() == 0) // empty line
		{
			++lineno;
			continue;
		}

		// the first token is the command name, pop it off and then process the
		// rest of the tokens
		string command { molstat::to_lower(tokens.front()) };
		tokens.pop();

		// go through the supported commands
		if(command == "endmodel")
		{
			// we're done here
			return ret;
		}
		else if(command == "model")
		{
			// the total lineno needs to be passed to the (sub)model reader
			const size_t mylineno{ lineno };

			// enter the (sub)model reader to process input lines until the
			// appropriate "endmodel" command is found
			ModelInformation model = readModel(input, output, ++lineno);

			// make sure there is a name (model type) specified
			if(tokens.size() == 0)
			{
				printError(output, mylineno, "No submodel type specified.");
			}
			else
			{
				// store the name of the model type
				model.name = molstat::to_lower(tokens.front());

				// set this to the parser's top model
				ret.submodels.emplace_back( move(model) );
			}
		}
		else if(command == "distribution")
		{
			// make sure there are tokens, if so, push the tokens
			if(tokens.size() < 2)
			{
				printError(output, lineno,
					"No distribution name and/or type specified.");
			}
			else
			{
				const string name{ tokens.front() };
				tokens.pop();

				// construct the random number distribution
				try
				{
					shared_ptr<const molstat::RandomDistribution> dist
						{ molstat::RandomDistributionFactory(move(tokens)) };
					ret.dists.emplace(name, dist);
				}
				catch(const invalid_argument &e)
				{
					// indent the error message
					printError(output, lineno,
						molstat::find_replace(e.what(), "\n", "\n   "));
				}
			}
		}
		else
		{
			printError(output, lineno,
				"Unknown model command: \"" + command + "\".");
		}

		// move to the next line
		++lineno;
	}

	// if we're here, we hit EOF before finding the appropriate endmodel command
	throw runtime_error("Missing \"endmodel\" command.");
}

std::shared_ptr<molstat::SimulateModel> SimulatorInputParse::constructModel(
	std::ostream &output,
	const std::map<std::string,
	               molstat::SimulateModelFactoryFunction> &models,
	ModelInformation &info)
{
	// see if the name specified is valid
	if(models.count(info.name) == 0)
		throw runtime_error("Unknown model: \"" + info.name + "\".");

	molstat::SimulateModelFactory factory{ models.at(info.name)() };

	// set the distributions, removing any distributions that aren't used
	{
		auto dist_iter = info.dists.cbegin();
		while(dist_iter != info.dists.cend())
		{
			bool used;
			factory.setDistribution(dist_iter->first, dist_iter->second, &used);
			if(!used)
			{
				// this distribution wasn't used... remove it from the list
				auto here = dist_iter;
				++dist_iter;
				info.dists.erase(here);
			}
			else
				++dist_iter;
		}
	}

	// add any submodels and remove any that aren't compatible/usable
	{
		auto submodel_iter = info.submodels.begin();
		while(submodel_iter != info.submodels.end())
		{
			try
			{
				// create the submodel
				shared_ptr<molstat::SimulateModel> submodel 
					{ constructModel(output, models, *submodel_iter) };

				// add the submodel
				factory.addSubmodel(submodel);

				// advance to the next submodel
				++submodel_iter;
			}
			catch(const exception &e)
			{
				output << "Error: " << e.what() << endl;

				// this submodel was unusable... remove it from the information
				auto here = submodel_iter;
				++submodel_iter;
				info.submodels.erase(here);
			}
		}
	}

	// get the model
	shared_ptr<molstat::SimulateModel> model{ nullptr };
	try
	{
		model = factory.getModel();
	}
	catch(const logic_error &e)
	{
		throw runtime_error(string(e.what()) + " When constructing model\n   " + 
			molstat::find_replace(info.to_string(), "\n", "\n   "));
	}

	return model;
}

std::size_t SimulatorInputParse::numTrials() const noexcept
{
	return trials;
}

std::string SimulatorInputParse::ModelInformation::to_string() const
{
	// first put in the name
	string ret{ name };
	ret += "\n   " + std::to_string(dists.size()) + " Distribution";
	if(dists.size() != 1)
		ret += 's';
	if(dists.size() > 0)
		ret += ':';

	// load in the distributions
	for(auto dist : dists)
	{
		ret += "\n      " + dist.first + " -> " + dist.second->info();
	}

	// submodel information
	for(auto submodel : submodels)
	{
		ret += "\n   Submodel type: ";

		// get the submodel info; indent it
		ret += molstat::find_replace(submodel.to_string(), "\n", "\n   ");
	}

	return ret;
}

void SimulatorInputParse::printState(std::ostream &output) const
{
	output << "Model type: " << top_model.to_string() << "\n\n";

	output << "Observables:\n";
	{
		for(auto obs_bin : obs_bins)
		{
			output << obs_bin.first << " -> " << obs_bin.second.first <<
				" (" << obs_bin.second.second->info() << ")\n";
		}
	}
	output << '\n';

	output << trials << " data point";
	if(trials != 1)
		output << 's';
	output << " will be simulated.\n";

	output << "Histogram Output File: " << histfilename << '\n';
}

std::string SimulatorInputParse::outputFileName() const
{
	return histfilename;
}

std::vector<std::shared_ptr<molstat::BinStyle>>
	SimulatorInputParse::getBinStyles() const
{
	vector<shared_ptr<molstat::BinStyle>> ret(obs_bins.size());

	for(auto obs_bin : obs_bins)
	{
		ret[obs_bin.first] = obs_bin.second.second;
	}

	return ret;
}