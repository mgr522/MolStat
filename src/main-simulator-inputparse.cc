/* This file is a part of MolStat, which is distributed under the Creative
   Commons Attribution-NonCommercial 4.0 International Public License.
   MolStat (c) 2014, Northwestern University. */
/**
 * \file main-simulator-inputparse.cc
 * \brief Definition of the SimulatorInputParse class for parsing the simulator
 *    input deck.
 *
 * \author Matthew G.\ Reuter
 * \date October 2014
 */

#include "main-simulator.h"
#include <iomanip>
#include <general/simulator_tools/simulator_exceptions.h>
#include <general/random_distributions/rng.h>
#include <electron_transport/simulator_models/transport_simulate_module.h>

using namespace std;

inline void SimulatorInputParse::printError(std::ostream &output,
	std::size_t lineno, std::string message)
{
	output << "Error on line " << setw(2) << lineno << ": " << message << endl;
}

void SimulatorInputParse::readInput(std::istream &input)
{
	std::size_t lineno{ 1 }; // line number

	// process the input deck
	// input evaluates to "false" when we hit EOF
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
		if(command == "bin" || command == "bin_x" || command == "bin_y")
		{
			// make sure we have some tokens for the binning style
			if(tokens.size() == 0)
			{
				printError(cout, lineno, "No binning information specified.");
			}
			else
			{
				// store the tokens for later
				if(command == "bin_y")
					bin_styles[1] = tokens;
				else
					bin_styles[0] = tokens;
			}
		}
		else if(command == "model")
		{
			// the total lineno needs to be passed to the model reader
			const size_t mylineno{ lineno };

			// enter the model reader to process input lines until the
			// appropriate "endmodel" command is found
			ModelInformation model = readModel(input, ++lineno);

			// make sure there is a name (model type) specified
			if(tokens.size() == 0)
			{
				printError(cout, mylineno, "No model type specified.");
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
			// make sure we have a name for the observable
			if(tokens.size() == 0)
			{
				printError(cout, lineno, "No observable specified.");
			}
			else {
				// store the name of the observable for later
				if(command == "observable_y")
					used_observables[1] = molstat::to_lower(tokens.front());
				else
					used_observables[0] = molstat::to_lower(tokens.front());
			}
		}
		else if(command == "output")
		{
			if(tokens.size() == 0)
			{
				printError(cout, lineno, "No output file name specified.");
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
				printError(cout, lineno, "Number of trials not specified.");
			}
			else
			{
				try
				{
					trials = molstat::cast_string<size_t>(tokens.front());
				}
				catch(const bad_cast &e)
				{
					printError(cout, lineno, "Unable to convert \"" + tokens.front() +
						"\" to a non-negative number.");
				}
			}
		}
		else
		{
			printError(cout, lineno, "Unknown command: \"" + command + "\".");
		}

		// move to the next line
		++lineno;
	}
}

SimulatorInputParse::ModelInformation SimulatorInputParse::readModel(
	std::istream &input, std::size_t &lineno)
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
			ModelInformation model = readModel(input, ++lineno);

			// make sure there is a name (model type) specified
			if(tokens.size() == 0)
			{
				printError(cout, mylineno, "No submodel type specified.");
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
				printError(cout, lineno,
					"No distribution name and/or type specified.");
			}
			else
			{
				const string name{ tokens.front() };
				tokens.pop();

				// construct the random number distribution
				shared_ptr<const molstat::RandomDistribution> dist{ nullptr };
				try
				{
					dist = molstat::RandomDistributionFactory(move(tokens));
					ret.dists.emplace(name, dist);
				}
				catch(const invalid_argument &e)
				{
					// indent the error message
					printError(cout, lineno,
						molstat::find_replace(e.what(), "\n", "\n   "));
				}
			}
		}
		else
		{
			printError(cout, lineno,
				"Unknown model command: \"" + command + "\".");
		}

		// move to the next line
		++lineno;
	}

	// if we're here, we hit EOF before finding the appropriate endmodel command
	throw runtime_error("Missing \"endmodel\" command.");
}

std::unique_ptr<molstat::Simulator> SimulatorInputParse::createSimulator()
{
	// the map of model instantiators
	map<string, molstat::SimulateModelFactoryFunction> models;

	// the map of observable indexes
	map<string, molstat::ObservableIndex> observables;

	// load model and observable names
	// the general syntax for loading models is
	// models.emplace( to_lower(name),
	//                 molstat::GetSimulateModelFactory<model_type>() );
	//
	// and likewise for observables,
	// observables.emplace( to_lower(name),
	//                      GetObservableIndex<observable_type>() );
	molstat::transport::load_models(models);
	molstat::transport::load_observables(observables);

	// make the model
	// if there are exceptions, let them pass up to the caller
	shared_ptr<molstat::SimulateModel> model
		{ constructModel(models, top_model) };
	
	// make the simulator
	unique_ptr<molstat::Simulator> sim{ new molstat::Simulator(model) };

	// set the observables
	for(const auto obs : used_observables)
	{
		try
		{
			auto obsindex = observables.at(obs.second);
			sim->setObservable(obs.first, obsindex);
		}
		catch(const out_of_range &e)
		{
			cout << "Unknown observable: \"" << obs.second << "\"." << endl;
		}
		catch(const logic_error &e)
		{
			cout << "Error setting observable " << obs.first << ".\n   " <<
				e.what() << endl;
		}
	}

	return sim;
}

std::shared_ptr<molstat::SimulateModel> SimulatorInputParse::constructModel(
	const std::map<std::string,
	               molstat::SimulateModelFactoryFunction> &models,
	const ModelInformation &info)
{
	// see if the name specified is valid
	if(models.count(info.name) == 0)
		throw runtime_error("Unknown model: \"" + info.name + "\".");

	molstat::SimulateModelFactory factory{ models.at(info.name)() };

	// set the distributions
	for(auto dist : info.dists)
	{
		factory.setDistribution(dist.first, dist.second);
	}

	// add any submodels
	for(auto submodelinfo : info.submodels)
	{
		shared_ptr<molstat::SimulateModel> submodel{ nullptr };

		// create the submodel
		try
		{
			submodel = constructModel(models, submodelinfo);

			// add the submodel
			try
			{
				factory.addSubmodel(submodel);
			}
			catch(const exception &e)
			{
				cout << "Error: " << e.what() << endl;
			}
		}
		catch(const exception &e)
		{
			cout << "Error: " << e.what() << endl;
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

std::string SimulatorInputParse::ModelInformation::to_string() const
{
	// first put in the name
	string ret{ name };
	ret += "\n   Distributions:";

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
	output << "--------------------------------------------------\n" <<
		"State of input parser.\n\n";

	output << "Model type: " << top_model.to_string() << "\n\n";

	output << "Observables:\n";
	{
		auto bin_iter = bin_styles.cbegin();
		auto obs_iter = used_observables.cbegin();
		while(bin_iter != bin_styles.cend() ||
			obs_iter != used_observables.cend())
		{
			if(bin_iter == bin_styles.cend() || obs_iter->first < bin_iter->first)
			{
				// observable specified but not binning style
				output << obs_iter->first << " -> " << obs_iter->second <<
					" (no binning information)\n";
				++obs_iter;
			}
			else if(obs_iter == used_observables.cend() ||
				bin_iter->first < obs_iter->first)
			{
				// bin specified but not observable
				output << bin_iter->first << " -> <missing observable> (" <<
					bin_iter->second.front() << ")\n";
				++bin_iter;
			}
			else /*if(obs_iter->first == bin_iter->first)*/
			{
				// same index
				output << obs_iter->first << " -> " << obs_iter->second <<
					" (" << bin_iter->second.front() << ")\n";
				++obs_iter;
				++bin_iter;
			}
		}
	}
	output << '\n';

	output << trials << " data point";
	if(trials != 1)
		output << 's';
	output << " will be simulated.\n";

	output << "Histogram Output File: " << histfilename << '\n';

	output << "--------------------------------------------------" << endl;
}